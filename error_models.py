import numpy as np
from numpy import sin, cos, exp
from numpy import sqrt
import global_file
from global_file import params
from functools import reduce
from qutip import tensor, qeye, ket2dm, basis, Qobj
from qutip.qip.gates import hadamard_transform, cnot
from itertools import product
from operator import mul
import pdb

def trace_out_middle_systems(state):
    reshaped_state = (state.full()).reshape([9, 81, 9, 9, 81, 9])
    reduced_state = np.einsum('ijkljn->ikln', reshaped_state)
    return Qobj(reduced_state.reshape([81, 81]))


def ptrace(state, system):
    d = (state.dims)[0][0]
    reshaped_state = (state.full()).reshape([d, d, d, d])
    if system == 0:
        reduced_state = np.einsum('jiki->jk', reshaped_state)
    elif system == 1:
        reduced_state = np.einsum('ijik->jk', reshaped_state)
    state = Qobj(reduced_state.reshape([d, d]))
    state.dims = [[d], [d]]
    return state


large_cnot = tensor(qeye(2), cnot(), qeye(2))
large_cnot.dims = [[2, 2, 2, 2], [2, 2, 2, 2]]

large_H = tensor(qeye(2), hadamard_transform(1), qeye(2), qeye(2))

M00 = tensor(qeye(2), ket2dm(basis(2, 0)), ket2dm(basis(2, 0)), qeye(2))
M01 = tensor(qeye(2), ket2dm(basis(2, 0)), ket2dm(basis(2, 1)), qeye(2))
M10 = tensor(qeye(2), ket2dm(basis(2, 1)), ket2dm(basis(2, 0)), qeye(2))
M11 = tensor(qeye(2), ket2dm(basis(2, 1)), ket2dm(basis(2, 1)), qeye(2))


# used to generate povms for four_qutrit_modes
MES1 = ket2dm((basis(81, 36) + basis(81, 28) + basis(81, 12) + basis(81, 4))/2)
MES2 = ket2dm((basis(81, 36) - basis(81, 28) + basis(81, 12) - basis(81, 4))/2)
MES3 = ket2dm((basis(81, 36) + basis(81, 28) - basis(81, 12) - basis(81, 4))/2)
MES4 = ket2dm((basis(81, 36) - basis(81, 28) - basis(81, 12) + basis(81, 4))/2)

# four_qutrit_modes povms
MQ1 = tensor(qeye(9), MES1, qeye(9))
MQ2 = tensor(qeye(9), MES2, qeye(9))
MQ3 = tensor(qeye(9), MES3, qeye(9))
MQ4 = tensor(qeye(9), MES4, qeye(9))

MQ1.dims = [[9, 9, 9, 9], [9, 9, 9, 9]]
MQ2.dims = [[9, 9, 9, 9], [9, 9, 9, 9]]
MQ3.dims = [[9, 9, 9, 9], [9, 9, 9, 9]]
MQ4.dims = [[9, 9, 9, 9], [9, 9, 9, 9]]


X = Qobj([[0, 1],
          [1, 0]])
Y = Qobj([[0, -1j],
          [1j, 0]])
Z = Qobj([[1, 0],
          [0, -1]])

large_X = tensor(X, qeye(2))
large_Y = tensor(Y, qeye(2))
large_Z = tensor(Z, qeye(2))

large_Z2 = tensor(qeye(2), Z)

Psi_plus = ket2dm(basis(4, 1) + basis(4, 2)).unit()
Psi_min = ket2dm(basis(4, 1) - basis(4, 2)).unit()

up_up = ket2dm(basis(4, 3))
mms = qeye(4)/4
mms.dims = [[2, 2], [2, 2]]


def transmissivity(distance):
    return np.exp(-distance / global_file.params.Latt)


def perform_BSM(state_A, state_B, node):
    if global_file.params.simulation_type == "bell_diagonal":
        b1A, b2A, b3A, b4A = state_A
        b1B, b2B, b3B, b4B = state_B

        b1 = b1A * b1B + b2A * b2B + b3A * b3B + b4A * b4B
        b2 = b1A * b2B + b2A * b1B + b3A * b4B + b4A * b3B
        b3 = b1A * b3B + b2A * b4B + b3A * b1B + b4A * b2B
        b4 = b1A * b4B + b2A * b3B + b3A * b2B + b4A * b1B

        state = [b1, b2, b3, b4]
        return state

    elif global_file.params.simulation_type == "IP":
        # Bell state measurements are always performed from qubit 2 in state1 to qubit 1 in state2

        state = large_cnot * tensor(state_A, state_B) * large_cnot
        state = depolarise_for_swap(state, node['CNOT_depol'])
        state = large_H * state * large_H

        state00 = (M00 * state * M00).ptrace([0, 3])
        state01 = (M01 * state * M01).ptrace([0, 3])
        state10 = (M10 * state * M10).ptrace([0, 3])
        state11 = (M11 * state * M11).ptrace([0, 3])

        state01 = large_X*state01*large_X
        state10 = large_Z*state10*large_Z
        state11 = large_Y*state11*large_Y

        state = state00 + state01 + state10 + state11
        return state

    elif global_file.params.simulation_type == "MP":
        state = tensor(state_A, state_B)
        state00 = MQ1 * state * MQ1
        state00 = trace_out_middle_systems(state00)
        state = state00/state00.tr()
        state.dims = [[9, 9], [9, 9]]
        return state


def perform_perfect_BSM(state_A, state_B):
    if global_file.params.simulation_type == "bell_diagonal":
        b1A, b2A, b3A, b4A = state_A
        b1B, b2B, b3B, b4B = state_B

        b1 = b1A * b1B + b2A * b2B + b3A * b3B + b4A * b4B
        b2 = b1A * b2B + b2A * b1B + b3A * b4B + b4A * b3B
        b3 = b1A * b3B + b2A * b4B + b3A * b1B + b4A * b2B
        b4 = b1A * b4B + b2A * b3B + b3A * b2B + b4A * b1B

        state = [b1, b2, b3, b4]
        return state

    elif global_file.params.simulation_type == "IP":
        # Bell state measurements are always performed from qubit 2 in state1 to qubit 1 in state2

        state = large_cnot * tensor(state_A, state_B) * large_cnot
        state = large_H * state * large_H

        state00 = (M00 * state * M00).ptrace([0, 3])
        state01 = (M01 * state * M01).ptrace([0, 3])
        state10 = (M10 * state * M10).ptrace([0, 3])
        state11 = (M11 * state * M11).ptrace([0, 3])

        state01 = large_X*state01*large_X
        state10 = large_Z*state10*large_Z
        state11 = large_Y*state11*large_Y

        state = state00 + state01 + state10 + state11
        return state


def concatenate_dephasing(dephasing_parameters):
    """
    Returns the total amount of dephasing for a list of dephasing parameters, where the
    list is in order of the noise application.
    """
    dephasing = 1
    for i in range(0, len(dephasing_parameters)):
        dephasing = dephasing * \
            dephasing_parameters[i] + \
            (1 - dephasing) * (1 - dephasing_parameters[i])
    return dephasing


def concatenate_depolarising(depolarising_parameters):
    """
    Returns the total amount of depolarising for a list of depolarising parameters.
    """
    return reduce(mul, depolarising_parameters)


def round_decoherence_parameters(x):
    if 1 < x < 1+(1e-6):
        return 1
    elif 1+1e-6 < x:
        raise ValueError("Gamma, tau or delta are greater than unity!")
    else:
        return x


def decohere_state(state, gamma, tau, delta = 1, system=1):
    gamma = round_decoherence_parameters(gamma)
    tau = round_decoherence_parameters(tau)
    delta = round_decoherence_parameters(delta)
    state = depolarise_state(state, gamma, system)
    state = dephase_state(state, tau, system)
    if global_file.params.simulation_type == "IP":
        if delta < 1:  # Don't apply damping if delta ==1
            state = damp_state(state, delta, system)
    return state


def depolarise_state(state, gamma, system=1):
    if global_file.params.simulation_type == "bell_diagonal":
        return [1 / 4 + (b - 1 / 4) * gamma for b in state]
    elif global_file.params.simulation_type == "IP":
        if system == 1:
            return (1-gamma) * tensor(qeye(2)/2, state.ptrace(1)) + gamma*state
        elif system == 2:
            return (1-gamma) * tensor(state.ptrace(0), qeye(2)/2) + gamma*state
        else:
            raise ValueError("No system specified to apply noise to")
    elif global_file.params.simulation_type == "MP":
        if system == 1:
            return (1-gamma) * tensor(qeye(9)/9, state.ptrace(1)) + gamma*state
        elif system == 2:
            return (1-gamma) * tensor(state.ptrace(0), qeye(9)/9) + gamma*state
        else:
            raise ValueError("No system specified to apply noise to")


def depolarise_for_swap(state, gamma):
    if global_file.params.simulation_type == "IP":
        return (1-gamma) * tensor(state.ptrace(0), mms, state.ptrace(3)) + gamma*state
    elif global_file.params.simulation_type == "MP":
        return state
    else:
        raise ValueError("Not implemented for bell diagonal formalism.")


def depolarise_for_distillation(state, gamma, node):
    if global_file.params.simulation_type == "IP":
        state2 = state
        if node == 1:
            state3 = (1-gamma) * tensor(mms, state2.ptrace([2, 3])) + gamma * state2
        elif node == 2:
            state3 = (1 - gamma) * tensor(state2.ptrace([0, 1]), mms) + gamma * state2
        else:
            raise ValueError("No/invalid node specified")
        return state3.permute([0, 2, 1, 3])
    elif global_file.params.simulation_type == "MP":
        return state
    else:
        raise ValueError("Not implemented for bell diagonal formalism.")


def dephase_state(state, tau, system=1):
    if global_file.params.simulation_type == "bell_diagonal":
        b1, b2, b3, b4 = state
        return [b1 * tau + (1 - tau) * b4, b2 * tau + (1 - tau) *
            b3, b3 * tau + (1 - tau) * b2, b4 * tau + (1 - tau) * b1]
    elif global_file.params.simulation_type == "IP":
        if system == 1:
            return tau * state + (1-tau) * (large_Z*state*large_Z)
        elif system == 2:
            return tau * state + (1-tau) * (large_Z2*state*large_Z2)
        else:
            raise ValueError("No system specified to apply noise to")
    elif global_file.params.simulation_type == "MP":
        no_coherence = Qobj(np.diag(state.diag()))
        no_coherence.dims = [[9, 9], [9, 9]]
        state = tau*state + (1-tau)*no_coherence
        return state


def damp_state(state, delta, system=1):
    if global_file.params.simulation_type == "IP":
        E_0 = Qobj([[1, 0],
                    [0, np.sqrt(delta)]])
        E_1 = Qobj([[0, np.sqrt(1-delta)],
                    [0, 0]])
        if system == 1:
            E_0 = tensor(E_0, qeye(2))
            E_1 = tensor(E_1, qeye(2))
        elif system == 2:
            E_0 = tensor(qeye(2), E_0)
            E_1 = tensor(qeye(2), E_1)
        else:
            raise ValueError("Qubit specified is not valid.")
        return E_0*state*E_0.dag() + E_1*state*E_1.dag()
    elif global_file.params.simulation_type == "MP":
        E_0 = Qobj([[1, 0, 0],
                    [0, np.sqrt(delta), 0],
                    [0, 0, delta]])
        E_1 = Qobj([[0, np.sqrt(1-delta), 0],
                    [0, 0, np.sqrt(2*delta*(1-delta))],
                    [0, 0, 0]])
        E_2 = Qobj([[0, 0, 1-delta],
                    [0, 0, 0],
                    [0, 0, 0]])
        if system == 1:
            A_0 = tensor(E_0, E_0, qeye(3))
            A_1 = tensor(E_0, E_1, qeye(3))
            A_2 = tensor(E_0, E_2, qeye(3))
            A_3 = tensor(E_1, E_0, qeye(3))
            A_4 = tensor(E_1, E_1, qeye(3))
            A_5 = tensor(E_1, E_2, qeye(3))
            A_6 = tensor(E_2, E_0, qeye(3))
            A_7 = tensor(E_2, E_1, qeye(3))
            A_8 = tensor(E_2, E_2, qeye(3))
        elif system == 2:
            A_0 = tensor(qeye(3), E_0, E_0)
            A_1 = tensor(qeye(3), E_0, E_1)
            A_2 = tensor(qeye(3), E_0, E_2)
            A_3 = tensor(qeye(3), E_1, E_0)
            A_4 = tensor(qeye(3), E_1, E_1)
            A_5 = tensor(qeye(3), E_1, E_2)
            A_6 = tensor(qeye(3), E_2, E_0)
            A_7 = tensor(qeye(3), E_2, E_1)
            A_8 = tensor(qeye(3), E_2, E_2)
        else:
            raise ValueError("System specified is not valid.")
        return A_0*state*A_0.dag() + A_1*state*A_1.dag() + A_2*state*A_2.dag() + \
               A_3*state*A_3.dag() + A_4*state*A_4.dag() + A_5*state*A_5.dag() + \
               A_6*state*A_6.dag() + A_7*state*A_7.dag() + A_8*state*A_8.dag()
    elif global_file.params.simulation_type == "bell_diagonal":
        print("Damping is not implemented for bell diagonal simulations.")
        return state


def average_exp(p, c, r):
    if p == 0:
        return 0
    if global_file.params.worst_case:
        return exp(-c * (r - 1))
    else:
        x = np.power(1-p, r)
        y = np.exp(c)
        return ((np.exp(c*(1-r))-y*x)*p) / ((1 + y*(-1 + p))*(1-x))


def apply_dark_count_noise(state, p, dcs=None, system=2):
    if dcs is None:
        dcs = global_file.params.dark_counts
    dc = dcs * global_file.params.integration_time
    alpha = p * np.exp(-dc) / (1 - (1 - p) * np.exp(-2 * dc))
    return depolarise_state(state, alpha, system=system)


def calc_rho_out_double_click(theta, eta, F_prep, opu_noise, Fm):
    a = opu_noise * (F_prep**2 + (1-F_prep)**2) + 2*F_prep * (1-F_prep) * (1-opu_noise)
    b = (1-opu_noise)*(F_prep**2 + (1-F_prep)**2) + 2*F_prep * (1-F_prep) * opu_noise
    rho_0 = ((2*sin(theta)**2)/(2-eta * cos(theta)**2)) * (a*Psi_plus + b * Psi_min) \
            + (((cos(theta)**2)*(2-eta))/(2-eta*cos(theta)**2)) * ket2dm(basis(4, 3))
    rho_2 = (1/(1-eta * cos(theta)**2)**2) * (sin(theta)**4 * ket2dm(basis(4, 0))
                                              + ((1-eta)*(cos(theta)**2) * sin(theta)**2) *
                                              (ket2dm(basis(4, 1)) + ket2dm(basis(4, 2)))
                                              + ((1-eta)**2)*cos(theta)**4 * ket2dm(basis(4, 3)))

    pd = global_file.params.dark_counts*global_file.params.integration_time
    p0 = (eta * cos(theta)**2)*(1-(eta/2)*cos(theta)**2)
    p2 = (1-eta * cos(theta)**2)**2

    Y = (2*p0) * (1-pd) + 2 * p2 * pd * (1-pd)

    rho_out = ((2*p0) * (1 - pd)*rho_0 + 2 * p2 * pd * (1-pd) * rho_2)/Y
    rho_out.dims = [[2, 2], [2, 2]]

    rho_AB = Fm**2*rho_out + (1-Fm)*Fm * (tensor(qeye(2)/2, rho_out.ptrace(1))
                                          + tensor(rho_out.ptrace(0), qeye(2)/2)) + ((1-Fm)**2) * mms
    rho_AB = large_X*rho_AB*large_X
    return rho_AB, Y


def create_initial_state_MP(eta1, eta2, eta3, eta4, Ns):
    p0 = 1/(Ns+1)**2
    p1 = 2*Ns/((Ns+1)**3)
    p2 = 1-p0-p1

    p_click = (2*eta3**2*p2 + eta2*eta3*(p1 + 2*p2)*(3*p1 + 6*p2 - 8*eta3*p2) +
               2*eta2**2*p2*(1 + 8*eta3**2*p2 - 4*eta3*(p1 + 2*p2)))/24.
    p_photons = (eta1*eta4*(-2*eta3**2*(-2 + eta4)*p2*(p1 - (-2 + eta1)*p2) +
                            eta2*eta3*(p1 - 2*(-2 + eta1)*p2)*(3*p1 + 2*(-3 + 4*eta3)*(-2 + eta4)*p2) +
                            2*(-2 + eta1)*eta2**2*p2*((-1 + 4*eta3)*p1 + (1 - 8*eta3 + 8*eta3**2)*(-2 + eta4)*p2)))\
                /(2*eta3**2*p2 + eta2*eta3*(p1 + 2*p2)*(3*p1 + 6*p2 - 8*eta3*p2) + 2*eta2**2*p2*
                  (1 + 8*eta3**2*p2 - 4*eta3*(p1 + 2*p2)))

    state = np.zeros((81, 81), dtype=float)

    state[10, 10] = (eta1 * eta4 * p2 * (-3 * eta3 ** 2 * (-1 + eta4) * (p1 - 2 * (-1 + eta1) * p2) + eta2 * eta3 * (
                -3 * (-2 + eta1 + 2 * eta3 + eta4 - 2 * eta3 * eta4) * p1 - 4 * (-1 + eta1) * (-5 + 8 * eta3) * (
                    -1 + eta4) * p2) + (-1 + eta1) * eta2 ** 2 * (
                                                     (-3 + 6 * eta3) * p1 + 2 * (3 - 16 * eta3 + 16 * eta3 ** 2) * (
                                                         -1 + eta4) * p2))) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[11, 11] = (-2 * (-1 + eta1) * eta1 * eta2 * (-1 + eta3) * (
                -2 * eta3 + eta2 * (-1 + 3 * eta3)) * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[12, 12] = (eta1 * eta4 * (
                -6 * eta3 ** 2 * (-1 + eta4) * p2 * (p1 - 2 * (-1 + eta1) * p2) + 2 * (-1 + eta1) * eta2 ** 2 * p2 * (
                    3 * (-1 + 6 * eta3) * p1 + 2 * (3 - 32 * eta3 + 32 * eta3 ** 2) * (
                        -1 + eta4) * p2) + eta2 * eta3 * (
                            9 * p1 ** 2 - 6 * (5 * eta1 + 5 * (-2 + eta4) - 6 * eta3 * (-1 + eta4)) * p1 * p2 - 8 * (
                                -1 + eta1) * (-13 + 16 * eta3) * (-1 + eta4) * p2 ** 2))) / (6. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[12, 28] = (eta1 * eta2 * eta3 * eta4 * (-3 * p1 + 8 * (-1 + eta1) * (-1 + eta2) * p2) * (
                -3 * p1 + 8 * (-1 + eta3) * (-1 + eta4) * p2)) / (6. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[13, 13] = -(eta1 * eta4 ** 2 * p2 * (
                2 * (-1 + eta1) * eta2 ** 2 * (1 - 8 * eta3 + 10 * eta3 ** 2) * p2 + 3 * eta2 * eta3 * (
                    -1 + 2 * eta3) * (p1 - 4 * (-1 + eta1) * p2) - 3 * eta3 ** 2 * (p1 - 2 * (-1 + eta1) * p2))) / (
                                3. * (2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                                                  1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[13, 21] = (sqrt(2) * eta1 * (-1 + eta2) * eta3 * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2)))

    state[13, 29] = -(sqrt(2) * eta1 * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 * (
                -3 * p1 + 8 * (-1 + eta1) * (-1 + eta2) * p2)) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[13, 37] = -((eta1 * (-eta3 + eta2 * (-1 + 2 * eta3)) * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[15, 15] = (-2 * eta1 * eta2 * (-1 + eta3) * eta4 ** 2 * p2 * (
                -((-1 + eta1) * eta2 * p2) + eta3 * (3 * p1 + (-1 + eta1) * (-10 + 11 * eta2) * p2))) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[15, 31] = -(sqrt(2) * eta1 * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 * (
                -3 * p1 + 8 * (-1 + eta1) * (-1 + eta2) * p2)) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[15, 39] = (sqrt(2) * eta1 * eta2 * (-1 + eta3) * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2)))

    state[19, 19] = (-2 * eta1 ** 2 * (-1 + eta2) * eta3 * (-eta3 + eta2 * (-2 + 3 * eta3)) * (
                -1 + eta4) * eta4 * p2 ** 2) / (3. * (2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                                                                  1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[21, 13] = (sqrt(2) * eta1 * (-1 + eta2) * eta3 * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2)))

    state[21, 21] = (-2 * eta1 ** 2 * (-1 + eta2) * eta3 * eta4 * p2 * (
                -(eta3 * (-1 + eta4) * p2) + eta2 * (3 * p1 + (-10 + 11 * eta3) * (-1 + eta4) * p2))) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[21, 37] = -(sqrt(2) * eta1 ** 2 * (-1 + eta2) * eta2 * eta3 * eta4 * p2 * (
                -3 * p1 + 8 * (-1 + eta3) * (-1 + eta4) * p2)) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[22, 22] = (2 * eta1 ** 2 * (-1 + eta2) * eta3 * (-eta3 + eta2 * (-2 + 3 * eta3)) * eta4 ** 2 * p2 ** 2) / (
                3. * (2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                                  1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[22, 38] = (4 * eta1 ** 2 * (-1 + eta2) * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[24, 24] = (8 * eta1 ** 2 * (-1 + eta2) * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[24, 40] = (4 * eta1 ** 2 * (-1 + eta2) * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[28, 12] = (eta1 * eta2 * eta3 * eta4 * (-3 * p1 + 8 * (-1 + eta1) * (-1 + eta2) * p2) * (
                -3 * p1 + 8 * (-1 + eta3) * (-1 + eta4) * p2)) / (6. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[28, 28] = (eta1 * eta4 * (
                -6 * eta3 ** 2 * (-1 + eta4) * p2 * (p1 - 2 * (-1 + eta1) * p2) + 2 * (-1 + eta1) * eta2 ** 2 * p2 * (
                    3 * (-1 + 6 * eta3) * p1 + 2 * (3 - 32 * eta3 + 32 * eta3 ** 2) * (
                        -1 + eta4) * p2) + eta2 * eta3 * (
                            9 * p1 ** 2 - 6 * (5 * eta1 + 5 * (-2 + eta4) - 6 * eta3 * (-1 + eta4)) * p1 * p2 - 8 * (
                                -1 + eta1) * (-13 + 16 * eta3) * (-1 + eta4) * p2 ** 2))) / (6. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[29, 13] = -(sqrt(2) * eta1 * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 * (
            -3 * p1 + 8 * (-1 + eta1) * (-1 + eta2) * p2)) / (3. * (
            2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
            3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[29, 29] = (-2 * eta1 * eta2 * (-1 + eta3) * eta4 ** 2 * p2 * (
                -((-1 + eta1) * eta2 * p2) + eta3 * (3 * p1 + (-1 + eta1) * (-10 + 11 * eta2) * p2))) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[29, 37] = (sqrt(2) * eta1 * eta2 * (-1 + eta3) * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2)))

    state[30, 30] = (eta1 * eta4 * p2 * (-3 * eta3 ** 2 * (-1 + eta4) * (p1 - 2 * (-1 + eta1) * p2) + eta2 * eta3 * (
                -3 * (-2 + eta1 + 2 * eta3 + eta4 - 2 * eta3 * eta4) * p1 - 4 * (-1 + eta1) * (-5 + 8 * eta3) * (
                    -1 + eta4) * p2) + (-1 + eta1) * eta2 ** 2 * (
                                                     (-3 + 6 * eta3) * p1 + 2 * (3 - 16 * eta3 + 16 * eta3 ** 2) * (
                                                         -1 + eta4) * p2))) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[31, 15] = -(sqrt(2) * eta1 * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 * (
                -3 * p1 + 8 * (-1 + eta1) * (-1 + eta2) * p2)) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[31, 31] = -(eta1 * eta4 ** 2 * p2 * (
                2 * (-1 + eta1) * eta2 ** 2 * (1 - 8 * eta3 + 10 * eta3 ** 2) * p2 + 3 * eta2 * eta3 * (
                    -1 + 2 * eta3) * (p1 - 4 * (-1 + eta1) * p2) - 3 * eta3 ** 2 * (p1 - 2 * (-1 + eta1) * p2))) / (
                                3. * (2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                                                  1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[31, 39] = -((eta1 * (-eta3 + eta2 * (-1 + 2 * eta3)) * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[31, 55] = (sqrt(2) * eta1 * (-1 + eta2) * eta3 * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2)))

    state[33, 33] = (-2 * (-1 + eta1) * eta1 * eta2 * (-1 + eta3) * (
                -2 * eta3 + eta2 * (-1 + 3 * eta3)) * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[37, 13] = -((eta1 * (-eta3 + eta2 * (-1 + 2 * eta3)) * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[37, 21] = -(sqrt(2) * eta1 ** 2 * (-1 + eta2) * eta2 * eta3 * eta4 * p2 * (
                -3 * p1 + 8 * (-1 + eta3) * (-1 + eta4) * p2)) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[37, 29] = (sqrt(2) * eta1 * eta2 * (-1 + eta3) * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2)))

    state[37, 37] = -(eta1 ** 2 * eta4 * p2 * (2 * eta3 ** 2 * (-1 + eta4) * p2 + eta2 * eta3 * (
                -3 * p1 - 4 * (-3 + 4 * eta3) * (-1 + eta4) * p2) + eta2 ** 2 * ((-3 + 6 * eta3) * p1 + 2 * (
                3 - 12 * eta3 + 10 * eta3 ** 2) * (-1 + eta4) * p2))) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[38, 22] = (4 * eta1 ** 2 * (-1 + eta2) * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[38, 38] = (2 * eta1 ** 2 * eta2 * (-1 + eta3) * (
                -2 * eta3 + eta2 * (-1 + 3 * eta3)) * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[39, 15] = (sqrt(2) * eta1 * eta2 * (-1 + eta3) * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2)))

    state[39, 31] = -((eta1 * (-eta3 + eta2 * (-1 + 2 * eta3)) * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[39, 39] = -(eta1 ** 2 * eta4 * p2 * (2 * eta3 ** 2 * (-1 + eta4) * p2 + eta2 * eta3 * (
                -3 * p1 - 4 * (-3 + 4 * eta3) * (-1 + eta4) * p2) + eta2 ** 2 * ((-3 + 6 * eta3) * p1 + 2 * (
                3 - 12 * eta3 + 10 * eta3 ** 2) * (-1 + eta4) * p2))) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[39, 55] = -(sqrt(2) * eta1 ** 2 * (-1 + eta2) * eta2 * eta3 * eta4 * p2 * (
                -3 * p1 + 8 * (-1 + eta3) * (-1 + eta4) * p2)) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[40, 24] = (4 * eta1 ** 2 * (-1 + eta2) * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[40, 40] = (2 * eta1 ** 2 * (eta2 + eta3 - 2 * eta2 * eta3) ** 2 * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[40, 56] = (4 * eta1 ** 2 * (-1 + eta2) * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[42, 42] = (2 * eta1 ** 2 * eta2 * (-1 + eta3) * (
                -2 * eta3 + eta2 * (-1 + 3 * eta3)) * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[42, 58] = (4 * eta1 ** 2 * (-1 + eta2) * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[55, 31] = (sqrt(2) * eta1 * (-1 + eta2) * eta3 * eta4 * sqrt(eta1 * eta2 * eta3 * eta4) * p1 * p2) / (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2)))

    state[55, 39] = -(sqrt(2) * eta1 ** 2 * (-1 + eta2) * eta2 * eta3 * eta4 * p2 * (
                -3 * p1 + 8 * (-1 + eta3) * (-1 + eta4) * p2)) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[55, 55] = (-2 * eta1 ** 2 * (-1 + eta2) * eta3 * eta4 * p2 * (
                -(eta3 * (-1 + eta4) * p2) + eta2 * (3 * p1 + (-10 + 11 * eta3) * (-1 + eta4) * p2))) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[56, 40] = (4 * eta1 ** 2 * (-1 + eta2) * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[56, 56] = (8 * eta1 ** 2 * (-1 + eta2) * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[57, 57] = (-2 * eta1 ** 2 * (-1 + eta2) * eta3 * (-eta3 + eta2 * (-2 + 3 * eta3)) * (
                -1 + eta4) * eta4 * p2 ** 2) / (3. * (2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                                                                  1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[58, 42] = (4 * eta1 ** 2 * (-1 + eta2) * eta2 * (-1 + eta3) * eta3 * eta4 ** 2 * p2 ** 2) / (3. * (
                2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                            1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state[58, 58] = (2 * eta1 ** 2 * (-1 + eta2) * eta3 * (-eta3 + eta2 * (-2 + 3 * eta3)) * eta4 ** 2 * p2 ** 2) / (
                3. * (2 * eta3 ** 2 * p2 + eta2 * eta3 * (p1 + 2 * p2) * (
                    3 * p1 + 6 * p2 - 8 * eta3 * p2) + 2 * eta2 ** 2 * p2 * (
                                  1 + 8 * eta3 ** 2 * p2 - 4 * eta3 * (p1 + 2 * p2))))

    state = Qobj(state / p_photons)
    state.dims = [[9, 9], [9, 9]]

    return state, 4 * p_click * p_photons


def project_from_MP_to_IP(state):
    """
    Projects state from a 81x81 matrix (living in a space spanned by (00, 01, 02, 10, 11, 12, 20, 21, 22)
    to a 4x4 matrix (which is spanned by (00, 01, 10, 11))
    :param state: 81x81 state to be projected
    :return state_new, p: resultant state and the probability that the projection measurement succeeds.
    """
    proj = ket2dm(basis(9, 1)) + ket2dm(basis(9, 3))
    proj = tensor(proj, proj)
    p = (proj*state).tr()
    state = ((proj*state*proj)/p).full()

    state_new = np.zeros([4, 4], dtype=complex)

    # tuples corresponding to the indices of the 4x4 subspace to be projected into
    tuples = [(30, 2), (12, 3), (28, 0), (10, 1)]
    for i in product(tuples, repeat=2):
        index1 = i[0][1]
        index2 = i[1][1]
        index3 = i[0][0]
        index4 = i[1][0]

        state_new[index1][index2] = state[index3][index4]
    state_new = Qobj(state_new)
    state_new.dims = [[2, 2], [2, 2]]
    return state_new, p
