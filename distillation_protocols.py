from error_models import *
from scheme_class import *
from network_functions import *
from qutip import tensor, qeye, ket2dm, basis, Qobj
from qutip.qip.gates import cnot

double_cnots = cnot(N=4, control=0, target=2) * cnot(N=4, control=1, target=3)
double_cnots.dims = [[2, 2, 2, 2], [2, 2, 2, 2]]

M00d = tensor(qeye(2), qeye(2), ket2dm(basis(2, 0)), ket2dm(basis(2, 0)))
M01d = tensor(qeye(2), qeye(2), ket2dm(basis(2, 0)), ket2dm(basis(2, 1)))
M10d = tensor(qeye(2), qeye(2), ket2dm(basis(2, 1)), ket2dm(basis(2, 0)))
M11d = tensor(qeye(2), qeye(2), ket2dm(basis(2, 1)), ket2dm(basis(2, 1)))

U1 = (1/np.sqrt(2))*Qobj([[1, -1j],
                          [-1j, 1]])
U2 = (1/np.sqrt(2))*Qobj([[1, 1j],
                          [1j, 1]])

U = tensor(U1, U2, U1, U2)
large_X = tensor(X, qeye(2))


def get_max_achieved_fidelity_and_its_time(path):
    schemes_in_path = path_dict[return_path_key(path)].values()  # retrieve all schemes in path_dict[path]
    if schemes_in_path:
        scheme_with_max_fidelity = max(schemes_in_path, key=lambda scheme: scheme.fidelity)
        return scheme_with_max_fidelity.fidelity, scheme_with_max_fidelity.time
    else:
        raise ValueError("No schemes found with a fidelity above fidelity threshold.")


def calc_distillation_threshold(path):
    x, _ = get_max_achieved_fidelity_and_its_time(path)
    # inverse of equation 3.4 in BSc thesis T. Schiet
    distillation_threshold = (-1 + 2 * x - 3 * np.sqrt(-1 + 6 * x - 4 * x * x)) / (8 * x - 10)
    return distillation_threshold


def add_distillation_protocols(G, path, subscheme1, subscheme2, distillation_threshold):
    """Iterates over all possible distillation protocols. Checks whether the
    two subschemes have close enough fidelity as a heuristic.
    """
    if subscheme1.fidelity > distillation_threshold and subscheme2.fidelity > distillation_threshold:
        if np.abs(subscheme1.fidelity - subscheme2.fidelity) < global_file.params.eps_dist:
            distillation_protocol(G, path, subscheme1, subscheme2)


def distillation_protocol(G, path, subscheme1, subscheme2):
    """Performs the DJEMPS distillation protocol between the two states stored in
    subscheme1 and subscheme2. Subscheme1 is the control, and subscheme2 is the target.
    """
    node1 = G.nodes[path[0]]
    node2 = G.nodes[path[-1]]

    T1_1, T2_1, Tamp_1 = node1['T1'], node1['T2'], node1['Tamp']
    T1_2, T2_2, Tamp_2 = node2['T1'], node2['T2'], node2['Tamp']

    # distance between the nodes at the end of path
    distance = get_distance_over_path(G, path)

    cc_t = distance / global_file.params.c  # Travel time for classical communication
    distillation_time = max(node1['distillation_time'], node2['distillation_time'])
    state1 = subscheme1.state

    # Check whether sequential distillation is used, and if so, apply the corresponding waiting time noise
    if global_file.params.seq_distillation:
        round_time = (subscheme1.time + subscheme2.time + cc_t)
        waiting_time = subscheme2.time
        state1 = decohere_state(state1, np.exp(-waiting_time / T1_1),
                                np.exp(-waiting_time / T2_1), np.exp(-waiting_time / Tamp_1), system=1)
        state1 = decohere_state(state1, np.exp(-waiting_time / T1_2),
                                np.exp(-waiting_time / T2_2), np.exp(-waiting_time / Tamp_2), system=2)

    else:
        round_time = max(subscheme1.time, subscheme2.time) + cc_t

    # min_fidelity is used as a threshold for the
    # fidelity of our new state - distillation should increase our
    # fidelity
    min_fidelity = min(subscheme1.fidelity, subscheme2.fidelity)

    if global_file.params.simulation_type == "bell_diagonal":
        b1A, b2A, b3A, b4A = state1
        b1B, b2B, b3B, b4B = subscheme2.state
        N = (b1A + b2A) * (b1B + b2B) + (b3A + b4A) * (b3B + b4B)

        b1 = b1A * b1B + b2A * b2B
        b2 = b3A * b4B + b4A * b3B
        b3 = b3A * b3B + b4A * b4B
        b4 = b2A * b1B + b1A * b2B
        state = [b1 / N, b2 / N, b3 / N, b4 / N]
    elif global_file.params.simulation_type == "IP":
        state2 = subscheme2.state
        state = U * tensor(state1, state2) * U.dag()
        state = double_cnots * state * double_cnots.dag()

        prob00 = (M00d * state).tr()

        prob11 = (M11d * state).tr()

        state00 = (M00d * state * M00d).ptrace([0, 1])

        state11 = (M11d * state * M11d).ptrace([0, 1])

        N = prob00 + prob11
        if np.imag(N) < 1e-8:
            N = np.real(N)
        else:
            raise ValueError("Probability has a too large complex contribution.")
        state = (state00 + state11) / N
    else:
        raise ValueError("No valid simulation type selected.")

    state = decohere_state(state, node1['CNOT_depol'], node1['CNOT_depol'], 1, system=1)
    state = decohere_state(state, node2['CNOT_depol'], node2['CNOT_depol'], 1, system=2)
    state = decohere_state(state, node1['meas_depol'], node1['meas_deph'], 1, system=1)
    state = decohere_state(state, node2['meas_depol'], node2['meas_deph'], 1, system=2)
    # Add two qubit gate noise

    cc_and_gate_time = cc_t + distillation_time

    state = decohere_state(state, np.exp(-cc_and_gate_time / T1_1),
                           np.exp(-cc_and_gate_time / T2_1), np.exp(-cc_and_gate_time / Tamp_1), system=1)
    state = decohere_state(state, np.exp(-cc_and_gate_time / T1_2),
                           np.exp(-cc_and_gate_time / T2_2), np.exp(-cc_and_gate_time / Tamp_2), system=2)

    p = N * subscheme1.prob * subscheme2.prob

    r_min = int(np.ceil(np.log(1 - global_file.params.min_global_prob) /
                        np.log(1 - p)))
    r_max = int(np.ceil(np.log(1 - global_file.params.max_global_prob) / np.log(1 - p)))

    # Iterate over the possible number of attempts
    for r in range(r_min, r_max + 1):
        time = r * round_time

        avg_dec_1_depol = average_exp(p, round_time / T1_1, r)
        avg_dec_1_deph = average_exp(p, round_time / T2_1, r)
        avg_dec_1_amp = average_exp(p, round_time / Tamp_1, r)

        avg_dec_2_depol = average_exp(p, round_time / T1_2, r)
        avg_dec_2_deph = average_exp(p, round_time / T2_2, r)
        avg_dec_2_amp = average_exp(p, round_time / Tamp_2, r)

        decohered_state = decohere_state(state, avg_dec_1_depol,
                                         avg_dec_1_deph, avg_dec_1_amp)
        decohered_state = decohere_state(decohered_state, avg_dec_2_depol,
                                         avg_dec_2_deph, avg_dec_2_amp)

        if calc_fidelity(decohered_state) <= min_fidelity:
            # stop attempting storing for multiple r if
            # the fidelity drops below the starting fidelity
            break

        prob = 1 - np.power(1 - p, r)
        scheme = Scheme(path, 'd' + str(r), subscheme1, subscheme2, decohered_state, time, prob)

        add_scheme(scheme)
