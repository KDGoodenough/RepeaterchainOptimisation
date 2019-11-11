import numpy as np
import math
import global_file
from global_file import params, path_dict
import collections
from qutip import Qobj
from collections import defaultdict


def calc_fidelity(state):
    if global_file.params.simulation_type == "bell_diagonal":
        return state[0]
    elif global_file.params.simulation_type in ["IP", "MP"]:
        return np.real((state * global_file.params.target_state).tr())


# Define special states
max_corr_state = Qobj([[1/2, 0, 0, 0],
                       [0,   0, 0, 0],
                       [0,   0, 0, 0],
                       [0,   0, 0, 1/2]])
max_corr_state.dims = [[2, 2], [2, 2]]

max_mix_state = Qobj([[1/4, 0, 0, 0],
                      [0,   1/4, 0, 0],
                      [0,   0, 1/4, 0],
                      [0,   0, 0, 1/4]])
max_mix_state.dims = [[2, 2], [2, 2]]


class Scheme:
    """
    Holds the stored schemes and relevent attributes.
    state - state holds either the bell coefficients or the density matrix
    norm_f - normalised fidelity over the length of path, used for heuristics
    norm_fs - returns the norm_f of subscheme1 and subscheme2
    """

    def __init__(self, path, protocoltype, subscheme1,
                 subscheme2, state, time, prob, supp_data=None):
        self.path = path
        self.protocoltype = protocoltype
        self.subscheme1 = subscheme1
        self.subscheme2 = subscheme2
        self.time = time
        self.state = state
        self.prob = prob

        self.fidelity = calc_fidelity(self.state)

        self.norm_f = self.fidelity**(1 / (len(self.path)))
        self.supp_data = supp_data

        if subscheme1:  # Checks if the protocol is over a non-elementary link
            self.subpath1, self.subpath2 = subscheme1.path, subscheme2.path
            # depth is a measure of how many subschemes were used for the creation of the scheme
            self.depth = max(self.subscheme1.depth, self.subscheme2.depth) + 1
            self.norm_fs = self.subscheme1.norm_f, self.subscheme2.norm_f
        else:
            self.subpath1, self.subpath2 = None, None
            self.depth = 1

    def print_parameters(self):
        print('F = ', self.fidelity, ', T = ', self.time)


def clear_schemes(path=None):
    global path_dict
    """
    Removes all stored schemes over path. If no argument is given or path
    equals None, remove all schemes.
    """
    # Currently removes all schemes by iterating over them. path_dict cannot be
    # set to an empty dictionary straightaway, because of some issues relating
    # to global variables
    if path is not None:
        path_dict[path] = Schemes()
    else:
        for key in list(path_dict.keys()):
            path_dict.pop(key)


class Schemes(collections.MutableMapping):
    """A dictionary that checks whether the key first exists,
        then adds it if it does not exist. If it exists,
        it checks whether the saved value has a lower time than
        the scheme about to be added. If so it does not add it."""

    def __init__(self, *args, **kwargs):
            """Use the object dict"""
            self.store = dict(*args, **kwargs)
            self.highest_f_schemes = {}  # tuples of fidelities and time
            self.schemes_with_fixed_prob = defaultdict(dict)

    def __setitem__(self, key, scheme):
        if scheme.fidelity > global_file.params.threshold_fidelity \
                and scheme.prob > global_file.params.min_global_prob:
            if key in self and self.schemes_with_fixed_prob[key[1]].get(key[0]):
                if self.store[key].time > scheme.time:
                    self.store[key] = scheme
                    self.schemes_with_fixed_prob[key[1]].update({key[0]: scheme})
            else:
                self.store[key] = scheme
                self.schemes_with_fixed_prob[key[1]].update({key[0]: scheme})

    def __getitem__(self, key):
        return self.store[key]

    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def add_scheme_to_Schemes(self, scheme):
        key = get_key(scheme, global_file.params.epsilon_f, global_file.params.epsilon_p)
        self[key] = scheme

    def prune_schemes(self):
        """
        # removes sub-optimal schemes. There are two pruning methods, total_pruning and no-total_pruning
        # No-total_pruning is used in the paper.
        """
        if global_file.params.total_pruning:
            # Compare every scheme stored with one another for sub-optimality. Runs in |schemes|^4.
            def is_suboptimal(scheme1, schemes):
                check_list = []

                key1 = get_key(scheme1, params.epsilon_f, params.epsilon_p)
                for scheme2 in schemes:
                    key2 = get_key(scheme1, params.epsilon_f, params.epsilon_p)
                    if key1[0] <= key2[0] and key1[1] <= key2[1] and scheme1.time > scheme2.time:
                        check_list.append(True)
                    else:
                        check_list.append(False)

                return check_list

            optimal_schemes = []
            all_schemes = [x[1] for x in self.store.items()]
            for _, scheme1 in self.store.items():
                if not any(is_suboptimal(scheme1, all_schemes)):
                    optimal_schemes.append(scheme1)

            self.store = {}
            for scheme in optimal_schemes:
                keys = get_key(scheme, params.epsilon_f, params.epsilon_p)
                self.store[keys] = scheme
        else:
            # Only consider pruning schemes within a fixed rounded success probability (i.e. for a given pslice)
            # this runs in |schemes|^2 time.
            self.store = {}
            schemes_with_fixed_prob2 = {}
            for key, pslice in sorted(self.schemes_with_fixed_prob.items()):  # iter schemes with fixed p tilde
                for index, key2 in enumerate(sorted(pslice, reverse=True)):
                    if index == 0:  # the scheme with highest fidelity should have highest time
                        pruned_dict = {key2: pslice[key2]}
                        previous_time = pslice[key2].time
                    else:
                        # The time should be monotonically decreasing, so if time is smaller, store it,
                        #  and set previous_time. Otherwise not.
                        if pslice[key2].time < previous_time:
                            pruned_dict[key2] = pslice[key2]
                            previous_time = pslice[key2].time
                schemes_with_fixed_prob2[key] = pruned_dict

            # Put all the pruned data back
            for key, pslice in schemes_with_fixed_prob2.items():
                for key2, scheme in pslice.items():
                    self.store[(key2, key)] = scheme


def get_key(scheme, epsilon_f, epsilon_p):
    """
    :param scheme: scheme to get key for
    :param epsilon_f: coarse-graining parameter for the fidelity
    :param epsilon_p: coarse-graining parameter for the probability
    :return tuple of coarse-grained fidelity and probability, which is used as a key:
    """
    fidelity_key = epsilon_round(scheme.fidelity, epsilon_f, global_file.params.threshold_fidelity)
    probability_key = epsilon_round(scheme.prob, epsilon_p, global_file.params.min_global_prob)
    return fidelity_key, probability_key


def epsilon_round(quantity, epsilon, lbound=0):
    """
    Rounds quantity down to the nearest multiple of epsilon if epsilon > 0.
    Returns quantity if epsilon = 0.
    """
    if epsilon > 0:
        return int(round((quantity - lbound) / epsilon - 1 / 2))
    elif epsilon == 0:
        return quantity
    else:
        raise ValueError('epsilon does not have a valid value.')


def add_scheme(scheme):
    """
    Adds scheme to stored data/path_dict. If uniform_sim, the key
    is the length of the path, otherwise the path is used as key.
    :param scheme: scheme to be added.
    :return:
    """
    if global_file.params.uniform_sim:
        key = len(scheme.path)
    else:
        key = scheme.path

    if key not in path_dict:
        path_dict[key] = Schemes()
    path_dict[key].add_scheme_to_Schemes(scheme)


def add_schemes(schemes):
    """Add a list of schemes. Mainly used for testing purposes."""
    for scheme in schemes:
        add_scheme(scheme)


def print_schemes_in_path(path):
    """Nicely prints out schemes stored over path and their fidelities and times. """
    print('Schemes stored over path', path, 'are:')
    for scheme in path_dict[path]:
        print('F = ', scheme.fidelity, ', T = ', scheme.time)
    print('------')


def prune_scheme_from_path(path):
    """
    Performs pruning over schemes stored over path. Also reports the changes in the number of schemes stored.
    :param path: path for the schemes to be pruned from.
    :return:
    """
    if global_file.params.uniform_sim:
        key = len(path)
    else:
        key = path
    num_schemes_prepruning = len(path_dict[key])
    path_dict[key].prune_schemes()
    num_schemes_postpruning = len(path_dict[key])
    num_removed_schemes = num_schemes_prepruning - num_schemes_postpruning
    print("Number of pruned schemes is {}, from {} to {}.".format(num_removed_schemes,
                                                                  num_schemes_prepruning,
                                                                  num_schemes_postpruning))
