from error_models import *
from scheme_class import *
from network_functions import get_ordered_nodes, get_distance_over_path


def add_bell_swaps(G, path, subscheme1, subscheme2):
    """
    Add bell swaps using subscheme1 and subscheme2.
    """
    n1 = len(subscheme1.path)
    n2 = len(subscheme2.path)
    F1 = subscheme1.fidelity
    F2 = subscheme2.fidelity

    # Use the banded swapping heuristic, i.e. that the ratio of F1^(1/n1) and F2^(1/n2)
    # should be close to unity, where the n's are the path lengths.

    if np.abs(np.log2(F1) / n1 - np.log2(F2) / n2) < global_file.params.eps_swap \
            and F1*F2+(1-F1)*(1-F2)/3 > global_file.params.threshold_fidelity:
        if global_file.params.simulation_type == "IP":
            scheme_list = swap_entanglement_protocol(G, path, subscheme1, subscheme2)
        elif global_file.params.simulation_type == "MP":
            scheme_list = swap_entanglement_protocol_MP(G, path, subscheme1, subscheme2)
        return scheme_list
    else:
        return []


def elementary_link_generation_MP(G, node1_name, node2_name):
    """
    :param G: nx.Graph - graph on which the elementary link generation takes place
    :param node1_name: str - id of node1 involved in elementary link generation
    :param node2_name: str - id of node2 involved in elementary link generation
    """
    # Assign relevant experimental parameters from the nodes
    node1, node2 = G.nodes[node1_name], G.nodes[node2_name]
    T1_1, T1_2 = node1['T1'], node2['T1']
    distance = G[node1_name][node2_name]['distance']
    papp = node1['pem']**4 * node1['pps']**4 * node2['pdet']**4
    number_of_fibres = G[node1_name][node2_name]['#fibres']

    p = transmissivity(distance/2) * papp  # prob of single link succeeding

    round_time = distance / global_file.params.c + max(node1['T_prep'], node2['T_prep'])
    # get min and max considered N_s, see paper for more information.
    max_Ns = 1.5*(-3 + sqrt(5 + (2*sqrt(3*global_file.params.threshold_fidelity
                                        + global_file.params.threshold_fidelity**2))/
                            global_file.params.threshold_fidelity))/2

    Ns_list = list(np.arange(0.0002, max_Ns, 0.0001))
    for Ns in Ns_list:
        state, p_succ = create_initial_state_MP(eta1=papp, eta2=p, eta3=p, eta4=papp, Ns=Ns)
        p_succ = 1 - np.power(1 - p_succ, number_of_fibres)

        # find the r's that achieve at least a probability of p_min and p_max
        if p_succ < 1/(50*number_of_fibres):
            r_min = 1
            r_max = 1
        elif p_succ != 1:
            r_min = int(np.ceil(np.log(1 - global_file.params.min_global_prob) /
                                np.log(1 - p_succ)))
            r_max = int(np.ceil(np.log(1 - global_file.params.max_global_prob) / np.log(1 - p_succ)))
        else:
            r_min = 1
            r_max = 1

        # calculate the range of r's to consider between r_min and r_max
        r_diff = r_max - r_min
        prob_diff = (global_file.params.max_global_prob - global_file.params.min_global_prob) / \
                    global_file.params.epsilon_p
        factor = 500

        if r_diff > factor * prob_diff:
            r_range = range(r_min, r_max + 1, int(r_diff/(prob_diff * factor)))
        else:
            r_range = range(r_min, r_max + 1)

        for r in r_range:
            # extra noise and probability due to repeating the protocol r times
            time = r * round_time
            avg_losses = average_exp(p_succ, round_time * (1/T1_1+1/T1_2), r)
            path = (node1_name, node2_name)
            prob = (1 - np.power(1 - p_succ, r)) * avg_losses  # effective success probability

            scheme = Scheme(path, 'Ns=' + str(Ns)+ ', r=' + str(r), None,
                            None, state, time, prob)  # assign None in scheme since there are no subschemes
            add_scheme(scheme)


def elementary_link_generation_MP_to_IP(G, node1_name, node2_name):
    """
        :param G: nx.Graph - graph on which the elementary link generation takes place
        :param node1_name: str - id of node1 involved in elementary link generation
        :param node2_name: str - id of node2 involved in elementary link generation
        """
    # Assign relevant experimental parameters from the nodes
    node1, node2 = G.nodes[node1_name], G.nodes[node2_name]
    T1_1, T2_1, Tamp_1 = node1['T1'], node1['T2'], node1['Tamp']
    T1_2, T2_2, Tamp_2 = node2['T1'], node2['T2'], node2['Tamp']
    distance = G[node1_name][node2_name]['distance']
    papp = node1['pem']**4 * node1['pps']**4 * node2['pdet']**4
    number_of_fibres = G[node1_name][node2_name]['#fibres']

    p = transmissivity(distance/2) * papp  # prob of single link succeeding
    # prob of at least single success

    round_time = distance / global_file.params.c + max(node1['T_prep'], node2['T_prep'])

    # get min and max considered N_s, see paper for more information.
    max_Ns = 1.5*(-3 + sqrt(5 + (2*sqrt(3*global_file.params.threshold_fidelity
                                        + global_file.params.threshold_fidelity**2))/
                            global_file.params.threshold_fidelity))/2.
    Ns_list = list(np.arange(0.00002, max_Ns, 0.0001))

    for Ns in Ns_list:
        state, p_succ = create_initial_state_MP(eta1=papp, eta2=p, eta3=p, eta4=papp, Ns=Ns)
        p_succ = 1 - np.power(1 - p_succ, number_of_fibres)
        state, p_proj = project_from_MP_to_IP(state)
        p_succ = p_succ * p_proj

        # find the r's that achieve at least a probability of p_min and p_max
        if 1 < p_succ < 1+1e-6:
            r_min = 1
            r_max = 1
        elif p_succ < 1/(50*number_of_fibres):
            r_min = 1
            r_max = 1
        elif p_succ != 1:
            r_min = int(np.ceil(np.log(1 - global_file.params.min_global_prob) /
                                np.log(1 - p_succ)))
            r_max = int(np.ceil(np.log(1 - global_file.params.max_global_prob) / np.log(1 - p_succ)))
        else:
            r_min = 1
            r_max = 1

        # calculate the range of r's to consider between r_min and r_max
        r_diff = r_max - r_min
        prob_diff = (global_file.params.max_global_prob - global_file.params.min_global_prob) \
                    / global_file.params.epsilon_p
        factor = 500

        if r_diff > factor * prob_diff:
            r_range = range(r_min, r_max + 1, int(r_diff/(prob_diff * factor)))
        else:
            r_range = range(r_min, r_max + 1)

        if p_succ > 1 / (20):
            for r in r_range:
                time = r * round_time
                # extra noise and probability due to repeating the protocol r times
                avg_dec_1_depol = average_exp(p_succ, round_time / T1_1, r)
                avg_dec_1_deph = average_exp(p_succ, round_time / T2_1, r)
                avg_dec_1_amp = average_exp(p_succ, round_time / Tamp_1, r)

                avg_dec_2_depol = average_exp(p_succ, round_time / T1_2, r)
                avg_dec_2_deph = average_exp(p_succ, round_time / T2_2, r)
                avg_dec_2_amp = average_exp(p_succ, round_time / Tamp_2, r)

                decohered_state = decohere_state(state, avg_dec_1_depol,
                                                 avg_dec_1_deph, avg_dec_1_amp, system=1)
                decohered_state = decohere_state(decohered_state, avg_dec_2_depol,
                                                 avg_dec_2_deph, avg_dec_2_amp, system=2)

                path = (node1_name, node2_name)
                prob = (1 - np.power(1 - p_succ, r))  # effective success probability

                scheme = Scheme(path, 'Ns=' + str(Ns)+ ', r=' + str(r), None, None,
                                decohered_state, time, prob)  # assign None in scheme since there are no subschemes
                add_scheme(scheme)


def single_click_protocol(G, node1_name, node2_name):
    """
        :param G: nx.Graph - graph on which the elementary link generation takes place
        :param node1_name: str - id of node1 involved in elementary link generation
        :param node2_name: str - id of node2 involved in elementary link generation
        """

    node1, node2 = G.nodes[node1_name], G.nodes[node2_name]
    T_prep = max(node1['T_prep'], node2['T_prep'])
    dcs = global_file.params.dark_counts
    T1_1, T2_1, Tamp_1 = node1['T1'], node1['T2'], node1['Tamp']
    T1_2, T2_2, Tamp_2 = node2['T1'], node2['T2'], node2['Tamp']

    distance = G[node1_name][node2_name]['distance']
    non_fibre_losses = (node1['pem'] * node2['pem']
                      * node1['pdet'] * node2['pdet']
                      * node1['pps'] * node2['pps'])

    eta = transmissivity(distance/2) * np.sqrt(non_fibre_losses)  # prob of single link succeeding
    # prob of at least single success

    round_time = distance / global_file.params.c + T_prep

    F_prep = node1['F_prep']
    opu_noise = node1['opu_noise']
    Fm = node1['Fm']

    theta_list = estimate_theta(F_prep, eta, opu_noise, Fm, round_time, T1_1, T2_1, Tamp_1, T1_2, T2_2, Tamp_2)

    for theta in theta_list:
        state, p = calc_rho_out_double_click(theta, eta, F_prep, opu_noise, Fm)

        state = decohere_state(state, np.exp(-round_time / T1_1),
                               np.exp(-round_time / T2_1), np.exp(-round_time / Tamp_1), system=1)
        state = decohere_state(state, np.exp(-round_time / T1_2),
                               np.exp(-round_time / T2_2), np.exp(-round_time / Tamp_2), system=2)
        state = apply_dark_count_noise(state, eta, dcs, system=2)

        # find the r's that achieve at least a probability of p_min and p_max
        if p > 0:
            r_min = int(
                np.ceil(
                    np.log(
                        1 -
                        global_file.params.min_global_prob) /
                    np.log(
                        1 -
                        p)))

            r_max = int(np.ceil(np.log(1 - global_file.params.max_global_prob) / np.log(1 - p)))
        else:
            r_min = 1
            r_max = 1

        # calculate the range of r's to consider between r_min and r_max
        r_range = range(r_min, r_max + 1, int(np.ceil((r_max-r_min)/200)))

        for r in r_range:
            # extra noise and probability due to repeating the protocol r times
            time = r * round_time
            avg_dec_1_depol = average_exp(p, round_time / T1_1, r)
            avg_dec_1_deph = average_exp(p, round_time / T2_1, r)
            avg_dec_1_amp = average_exp(p, round_time / Tamp_1, r)

            avg_dec_2_depol = average_exp(p, round_time / T1_2, r)
            avg_dec_2_deph = average_exp(p, round_time / T2_2, r)
            avg_dec_2_amp = average_exp(p, round_time / Tamp_2, r)

            decohered_state = decohere_state(state, avg_dec_1_depol,
                                             avg_dec_1_deph, avg_dec_1_amp, system=1)
            decohered_state = decohere_state(decohered_state, avg_dec_2_depol,
                                             avg_dec_2_deph, avg_dec_2_amp, system=2)
            path = (node1_name, node2_name)
            prob = 1 - np.power(1 - p, r)  # effective success probability
            scheme = Scheme(path, "sc {:.4f}".format(theta) + ' - ' + str(r), None,
                            None, decohered_state, time, prob)
            # assign None in scheme since there are no subschemes
            add_scheme(scheme)



def estimate_theta(F_prep, eta, opu_noise, Fm, round_time, T1_1, T2_1, Tamp_1, T1_2, T2_2, Tamp_2):
    """Estimates the values of theta to consider for the optimisation"""
    def calc_fidelity_for_theta(F_prep, eta, opu_noise, Fm, theta, round_time, T1_1, T2_1, Tamp_1, T1_2, T2_2, Tamp_2):
        state, p = calc_rho_out_double_click(theta, eta, F_prep, opu_noise, Fm)

        state = decohere_state(state, np.exp(-round_time / T1_1),
                               np.exp(-round_time / T2_1), np.exp(-round_time / Tamp_1), system=1)
        state = decohere_state(state, np.exp(-round_time / T1_2),
                               np.exp(-round_time / T2_2), np.exp(-round_time / Tamp_2), system=2)

        r = int(np.ceil(np.log(1 - global_file.params.max_global_prob) / np.log(1 - p)))

        avg_dec_1_depol = average_exp(p, round_time / T1_1, r)
        avg_dec_1_deph = average_exp(p, round_time / T2_1, r)
        avg_dec_1_amp = average_exp(p, round_time / Tamp_1, r)

        avg_dec_2_depol = average_exp(p, round_time / T1_2, r)
        avg_dec_2_deph = average_exp(p, round_time / T2_2, r)
        avg_dec_2_amp = average_exp(p, round_time / Tamp_2, r)

        decohered_state = decohere_state(state, avg_dec_1_depol,
                                         avg_dec_1_deph, avg_dec_1_amp, system=1)
        decohered_state = decohere_state(decohered_state, avg_dec_2_depol,
                                         avg_dec_2_deph, avg_dec_2_amp, system=2)
        fidelity1 = calc_fidelity(decohered_state)

        return fidelity1

    theta_list = np.linspace(1/2, np.pi/2-0.01, 300)
    fidelity_list = []
    for theta in theta_list:
        fid = calc_fidelity_for_theta(F_prep, eta, opu_noise, Fm, theta, round_time, T1_1, T2_1, Tamp_1, T1_2, T2_2, Tamp_2)
        fidelity_list.append(fid)

    theta_list = [theta_list[i] for i in range(len(fidelity_list)-1)
                  if fidelity_list[i]>global_file.params.threshold_fidelity and fidelity_list[i]<fidelity_list[i+1]]

    return theta_list


def double_click_protocol(G, node1_name, node2_name):
    node1, node2 = G.nodes[node1_name], G.nodes[node2_name]
    T_prep = max(node1['T_prep'], node2['T_prep'])
    Fp = min(node1['F_prep'], node2['F_prep'])
    T1_1, T2_1, Tamp_1 = node1['T1'], node1['T2'], node1['Tamp']
    T1_2, T2_2, Tamp_2 = node2['T1'], node2['T2'], node2['Tamp']

    distance = G[node1_name][node2_name]['distance']
    non_fibre_losses = (node1['pem'] * node2['pem']
                        * node1['pdet'] * node2['pdet']
                        * node1['pps'] * node2['pps'])

    eta = transmissivity(distance/2) * np.sqrt(non_fibre_losses)
    p = (eta**2)/4

    round_time = 2*distance / global_file.params.c + T_prep

    state = global_file.params.target_state
    state = decohere_state(state, 1, 1-4*Fp+12*Fp**2-16*Fp**3+8*Fp**4, 1)

    if np.log(1-p)<0:
        r_min = int(
            np.ceil(
                np.log(
                    1 -
                    global_file.params.min_global_prob) /
                np.log(
                    1 -
                    p)))

        r_max = int(np.ceil(np.log(1 - global_file.params.max_global_prob) / np.log(1 - p)))
    else:
        r_min = 1
        r_max = 1

    r_range = range(r_min, r_max + 1, int(np.ceil((r_max-r_min)/200)))

    for r in r_range:
        # extra noise and probability due to repeating the protocol r times
        time = r * round_time

        avg_dec_1_depol = average_exp(p, round_time / T1_1, r)
        avg_dec_1_deph = average_exp(p, round_time / T2_1, r)
        avg_dec_1_amp = average_exp(p, round_time / Tamp_1, r)

        avg_dec_2_depol = average_exp(p, round_time / T1_2, r)
        avg_dec_2_deph = average_exp(p, round_time / T2_2, r)
        avg_dec_2_amp = average_exp(p, round_time / Tamp_2, r)

        decohered_state = decohere_state(state, avg_dec_1_depol,
                                         avg_dec_1_deph, avg_dec_1_amp, system=1)
        decohered_state = decohere_state(decohered_state, avg_dec_2_depol,
                                         avg_dec_2_deph, avg_dec_2_amp, system=2)

        path = (node1_name, node2_name)
        prob = 1 - np.power(1 - p, r)  # effective success probability
        scheme = Scheme(path, "dc- " + str(r), None,
                        None, decohered_state, time, prob)
        # assign None in scheme since there are no subschemes
        add_scheme(scheme)



def swap_entanglement_protocol(G, path, subscheme1, subscheme2):
    """
    :param G: nx.Graph - graph on which the elementary link generation takes place
    :param path: tuple - tuple of ALL nodes involved in swap operations
    :param subscheme1: Scheme - scheme 1 used in the swapping protocol
    :param subscheme2: Scheme - scheme 2 used in the swapping protocol
            """

    path1, path2 = subscheme1.path, subscheme2.path

    if global_file.params.uniform_sim:
        # pick first node since the simulation is uniform anyway
        node1, node2, node3 = G.nodes[path[0]], G.nodes[path[0]], G.nodes[path[0]]
    else:
        node1, node2, node3 = get_ordered_nodes(G, path, path1, path2)

    T1_1, T2_1, Tamp_1 = node1['T1'], node1['T2'], node1['Tamp']
    T1_3, T2_3, Tamp_3 = node3['T1'], node3['T2'], node2['Tamp']
    state1, state2 = subscheme1.state, subscheme2.state
    t1, t2 = subscheme1.time, subscheme2.time

    distance1 = get_distance_over_path(G, subscheme1.path)
    distance2 = get_distance_over_path(G, subscheme2.path)
    cc_time = max([distance1, distance2]) / global_file.params.c

    round_time = max(t1, t2) + node2['bell_swap_time'] + cc_time

    state = perform_BSM(state1, state2, node2)
    state = decohere_state(state, node2['bsm_depol'], node2['bsm_deph'], 1)
    bsm_prob = node2['bsm_prob']
    p = subscheme1.prob * subscheme2.prob*bsm_prob

    # find the r's that achieve at least a probability of p_min and p_max
    r_min = int(
        np.ceil(
            np.log(1 - global_file.params.min_global_prob) /
            np.log(1 - p)))
    r_max = int(np.ceil(np.log(1 - global_file.params.max_global_prob) / np.log(1 - p)))

    scheme_list = []
    for r in range(r_min, r_max + 1):
        # extra noise and probability due to repeating the protocol r times
        avg_dec_1_depol = average_exp(p, round_time / T1_1, r)
        avg_dec_1_deph = average_exp(p, round_time / T2_1, r)
        avg_dec_1_amp = average_exp(p, round_time / Tamp_1, r)

        if global_file.params.uniform_sim:
            avg_dec_3_depol = avg_dec_1_depol
            avg_dec_3_deph = avg_dec_1_deph
            avg_dec_3_amp = avg_dec_1_amp
        else:
            avg_dec_3_depol = average_exp(p, round_time / T1_3, r)
            avg_dec_3_deph = average_exp(p, round_time / T2_3, r)
            avg_dec_3_amp = average_exp(p, round_time / Tamp_3, r)

        decohered_state = decohere_state(state, avg_dec_1_depol,
                                         avg_dec_1_deph, avg_dec_1_amp, system=1)
        decohered_state = decohere_state(decohered_state, avg_dec_3_depol,
                                         avg_dec_3_deph, avg_dec_3_amp, system=2)
        time = r * round_time
        prob = 1 - np.power(1 - p, r)

        scheme = Scheme(path, 'cSwap' + '-' + str(r), subscheme1, subscheme2,
                        decohered_state, time, prob)
        # stores scheme in scheme_list when multithreading,
        #  such that all the schemes can be combined afterwards
        if global_file.params.multi_threading:
            scheme_list.append(scheme)
        else:
            add_scheme(scheme)
    return scheme_list

def swap_entanglement_protocol_MP(G, path, subscheme1, subscheme2):
    """
        :param G: nx.Graph - graph on which the elementary link generation takes place
        :param path: tuple - tuple of ALL nodes involved in swap operations
        :param subscheme1: Scheme - scheme 1 used in the swapping protocol
        :param subscheme2: Scheme - scheme 2 used in the swapping protocol
                """
    path1, path2 = subscheme1.path, subscheme2.path

    if global_file.params.uniform_sim:
        # pick first node since the simulation is uniform anyway
        node1, node2, node3 = G.nodes[path[0]], G.nodes[path[0]], G.nodes[path[0]]
    else:
        node1, node2, node3 = get_ordered_nodes(G, path, path1, path2)

    T1_1, T1_2 = node1['T1'], node2['T1']
    #T1_2, T2_2 = node2['T1'], node2['T2']
    state1, state2 = subscheme1.state, subscheme2.state
    t1, t2 = subscheme1.time, subscheme2.time


    distance1 = get_distance_over_path(G, subscheme1.path)
    distance2 = get_distance_over_path(G, subscheme2.path)
    cc_time = max([distance1, distance2]) / global_file.params.c

    round_time = max(t1, t2) + node2['bell_swap_time'] + cc_time

    # Decoherence of state due to a single attempt
    state = perform_BSM(state1, state2, node2)
    state = decohere_state(state, node2['bsm_depol'], node2['bsm_deph'], 1)

    bsm_prob = node2['bsm_prob']
    p = subscheme1.prob * subscheme2.prob*bsm_prob
    # find the r's that achieve at least a probability of p_min and p_max
    if p<1:
        r_min = int(
            np.ceil(
                np.log(1 - global_file.params.min_global_prob) /
                np.log(1 - p)))
        r_max = int(np.ceil(np.log(1 - global_file.params.max_global_prob) / np.log(1 - p)))

    else:
        r_min = 1
        r_max = 1

    scheme_list = []
    for r in range(r_min, r_max + 1):
        # extra noise and probability due to repeating the protocol r times
        avg_losses_1 = average_exp(p, round_time / T1_1, r)
        avg_losses_2 = average_exp(p, round_time / T1_2, r)


        time = r * round_time
        prob = (1 - np.power(1 - p, r)) * avg_losses_1 * avg_losses_2

        scheme = Scheme(path, 'cSwap' + '-' + str(r), subscheme1, subscheme2,
                        state, time, prob)

        if global_file.params.multi_threading:
            scheme_list.append(scheme)
        else:
            add_scheme(scheme)
    return scheme_list
