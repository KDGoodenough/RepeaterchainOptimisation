"""
### Repeater optimisation algorithm ###
#######################################

Contains highest level functions for optimising repeater protocols.
All these functions use the optimise_schemes functions from connection_and_distillation.py.
"""

from plot_functions import *
from global_file import create_params_dict
from connection_and_distillation import *
import time
from itertools import product
import warnings
warnings.filterwarnings("error")
warnings.filterwarnings(
    "ignore", message="elementwise comparison failed; returning scalar instead, "
                      "but in the future will perform elementwise comparison")


def optimise_repeater_chain(distance, n, general_distillation,
                            plot=True, save=True, color=None, scatter=False,
                            label='Optimised protocols', symmetricoptimisation=False):
    """Perform optimisation over a repeater chain of a certain distance and n repeaters."""
    G = create_repeater_chain(distance, n)
    G = assign_parameters_to_nodes(G)

    optimise_schemes(G, general_distillation, symmetricoptimisation=symmetricoptimisation)

    if save:
        with open(label + '.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
                pickle.dump(path_dict, f)
    if plot:
        if scatter:
            plot_max_fidelity_with_time(G, 'A', 'B', color=color)
        else:
            plot_fidelity_vs_time_plot(G, 'A', 'B', label=label,
                                       color=color, connected=True)
    return G


def optimise_repeater_with_parameters_from_file(filename):
    global_file.params = create_params_dict('parameters/' + filename + '.txt')

    G = create_repeater_chain(global_file.params.distance, int(global_file.params.n))
    G = assign_parameters_to_nodes(G)

    optimise_schemes(G, global_file.params.general_distillation,
                     symmetricoptimisation=global_file.params.symmetricoptimisation)

    # with open('pkl_files/' + filename + '.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
    #     pickle.dump(path_dict, f)

    plot_fidelity_vs_time_plot(G, 'A', 'B', label=filename,
                               color=None, connected=True)
    return G


if __name__ == "__main__":
    filename = 'MP_params_set4_3_repeater_200km'
    optimise_repeater_with_parameters_from_file(filename=filename)

    plt.show()
