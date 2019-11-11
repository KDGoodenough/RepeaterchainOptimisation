""""
Contains all global constants needed for error models, and initialization of path_dict
"""

import qutip as qt
from numpy import pi
import numpy as np
from copy import copy
global path_dict


def bessel_quotient(v, x, eps):
    """
    Computes the quotient of I_{v+1}(x) / I_{v}(x).
    Uses algorithm described in 'Computation of Modified Bessel Functions and Their Ratios' D.E. Amos (1974)
    Used to calculate the noise due to phase uncertainty, see equation (A25) https://arxiv.org/pdf/1809.00364.pdf.

    Parameters
    ----------
    v : float
    x : float
    eps : float
        The precision
    """

    def _init_approx(v, x):
        # eq. (20a)
        approx = x / (v + (1 / 2) + np.sqrt((v + (3 / 2)) ** 2 + x ** 2))

        return approx

    def _recursion(v, x, r_num, r_denom):
        # eq. (20b)
        R = r_num / r_denom

        r_ret = x / (v + 1 + np.sqrt((v + 1) ** 2 + x ** 2 * R))

        return r_ret

    def _error(v, x, r_curr, r_prev):
        # eq. after Figure 1
        d = 1 / ((2 * (v + 1) / x) + r_prev)

        err = np.abs(d - r_curr) / r_curr

        return err

    if v < 0:
        raise ValueError("v needs to be non-negative")
    if x < 0:
        raise ValueError("x needs to be non-negative")
    if eps <= 0:
        raise ValueError("eps needs to be positive")
    if x == 0:
        return 0

    # Inital values
    err = eps + 1
    m = 1
    prev_diag = [_init_approx(v, x)]

    while err > eps:

        curr_diag = [_init_approx((v + m), x)]

        for i in range(m):
            curr_diag.append(_recursion((v + m - 1 - i), x, curr_diag[i], prev_diag[i]))

        err = _error(v, x, curr_diag[-1], curr_diag[-2])

        prev_diag = curr_diag

        m += 1

    return curr_diag[-1]


def convert_from_str(string):
    if string == "False":
        return False
    elif string == "True":
        return True
    elif string == "{}":
        return {}
    elif string == '"IP"' or string == '“IP"':
        return "IP"
    elif string == '"MP"' or string == '“MP"':
        return "MP"
    else:
        str_turned_to_float = float(string)
        if str_turned_to_float.is_integer():
            return int(str_turned_to_float)
        else:
            return str_turned_to_float


class Bunch(object):
    def __init__(self, adict):
        self.__dict__.update(adict)


params = dict()


def create_params_dict(filename):
    with open(filename) as f:
        filedata = [line[:line.rfind("#")].split() for line in f if not line.startswith("#") and line.strip()]

        params = {}
        for line in filedata:
            params[line[0]] = convert_from_str(line[2])
        params = Bunch(params)

        if params.simulation_type == "IP":
            target_state = qt.Qobj([[1 / 2, 0, 0, 1 / 2],
                                    [0, 0, 0, 0],
                                    [0, 0, 0, 0],
                                    [1 / 2, 0, 0, 1 / 2]])
            target_state.dims = [[2, 2], [2, 2]]
        elif params.simulation_type == "bell_diagonal":
            target_state = [1, 0, 0, 0]
        elif params.simulation_type == "MP":
            target_state = (qt.basis(81, 28) + qt.basis(81, 12))/(np.sqrt(2))
            target_state = qt.ket2dm(target_state)
            target_state.dims = [[9, 9], [9, 9]]

        params.target_state = target_state
        params.state_prep = copy(target_state)

        opu = (params.phase_uncertainty / 360) * 2 * pi  # optical phase uncertainty
        params.opu_noise = (bessel_quotient(0, 1 / opu ** 2, 1e-7) + 1) / 2  # equation (A25) https://arxiv.org/pdf/1809.00364.pdf

        if not hasattr(params, 'Fm'):
            params.Fm = params.F_g
        if not hasattr(params, 'meas_depol'):
            params.meas_depol = params.F_g
        if not hasattr(params, 'bsm_depol'):
            params.bsm_depol = params.F_g
        if not hasattr(params, 'CNOT_depol'):
            params.CNOT_depol = params.F_g
        return params

path_dict = {}
