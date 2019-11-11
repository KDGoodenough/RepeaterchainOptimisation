from global_file import path_dict
import global_file
import warnings
from connection_protocols import *
from distillation_protocols import *
from network_functions import *
from numpy import arange
import time

from multiprocessing import Pool

warnings.filterwarnings("ignore", message="elementwise comparison failed; returning scalar "
                                          "instead, but in the future will perform elementwise comparison")


def optimise_schemes(G, general_distillation=False, Alice='A', Bob='B', symmetricoptimisation=False):
    """

    :param G: nx.Graph. Graph over which to perform the optimisation. The nodes of G should have
        several attributes, pos for absolute position of the node, and error parameters.
    :param general_distillation: bool - whether to perform distillation over two of the
        same pair (=False) or any possible combination (=True)
    :param Alice: ID used for Alice
    :param Bob: ID used for Bob
    :param symmetricoptimisation: bool - indicates whether the bisection heuristic is to be applied or not
    :return:
    """

    clear_schemes()  # Removes schemes from possible previous optimisation
    if global_file.params.simulation_type not in ["bell_diagonal", "IP", "MP"]:
        raise Exception("No proper simulation type specified.")

    candidate_paths, _ = find_candidate_paths(G, Alice, Bob)

    G = assign_lengths(G)
    subpaths = find_all_subpaths(candidate_paths)
    max_length = get_max_length(candidate_paths)
    subpaths = init_elementary_links(G, subpaths, general_distillation)

    print('Finished elementary links.')
    k_range = create_k_range(symmetricoptimisation, max_length)

    for k, symmetric in k_range:
        subpaths = iterate_over_paths_of_length_k(
            G, subpaths, general_distillation, k, symmetric)
        print(np.round(k / max_length * 100, 2), "%")
    return candidate_paths


def create_k_range(symmetric, nodes):
    if symmetric:
        g, h = calc_g_h(nodes - 1)
        k_range1 = list(arange(3, g + 2))
        symmetric_range1 = [False]*len(k_range1)
        k_range2 = list((2**i)*g + 1 for i in range(1, h))
        k_range = k_range1 + k_range2
        symmetric_range2 = [True]*len(k_range2)
        symmetric_range = symmetric_range1 + symmetric_range2
        return iter(zip(k_range, symmetric_range))
    else:
        k_range = arange(3, nodes + 1)
        symmetric_range = [False]*len(list(k_range))
        return iter(zip(k_range, symmetric_range))


def calc_g_h(n):
    """
    :param n: int
    :return:
    returns integer g and h such that n = g^h and g is as small as possible."""

    d = 1
    while n % 2 == 0:
        n /= 2
        d += 1
        return int(n), int(d)


def iterate_over_paths_of_length_k(G, subpaths, general_distillation, k, symmetric):
    time1 = time.time()
    # get paths of length k from subpaths, and remove
    # them from subpaths for efficiency

    paths_of_length_k, subpaths = get_paths_of_length_k(subpaths, k)
    if global_file.params.uniform_sim:
        paths_of_length_k = [paths_of_length_k[0]]

    for path_of_length_k in paths_of_length_k:
        add_connection_protocols(G, path_of_length_k, symmetric)
        prune_scheme_from_path(path_of_length_k)
    time2 = time.time()
    print("Time for swap protocols = {0:10.5f}".format(time2 - time1))
    for path_of_length_k in paths_of_length_k:
        distillation_threshold = calc_distillation_threshold(path_of_length_k)
        for _ in range(global_file.params.distillation_rounds):
            for schemes in iterate_distillation_subschemes(path_of_length_k,
                                                           general_distillation):
                add_distillation_protocols(
                    G, path_of_length_k, schemes[0], schemes[1], distillation_threshold)
            prune_scheme_from_path(path_of_length_k)
    time3 = time.time()
    print("Time for dist protocols = {0:10.5f}".format(time3 - time2))
    print('--')
    return subpaths


def add_connection_protocols(G, path, symmetric):
    """
    Applies and stores the sufficiently good connection protocols of path.
    """
    if len(path) > 2:
        if global_file.params.multi_threading:
            data = [(G, path, bipartition, symmetric) for bipartition in get_bipartitions(path, symmetric)]
            with Pool(processes=None) as p:
                list_of_scheme_lists = p.map(sub_func, data)
                for list_of_schemes in list_of_scheme_lists:
                    add_schemes(list_of_schemes)
        else:
            for bipartition in get_bipartitions(path, symmetric):
                subpath1, subpath2 = bipartition[0], bipartition[1]
                for schemes in iterate_over_schemes_of_subpaths(subpath1, subpath2, symmetric):
                    add_bell_swaps(
                        G, path, schemes[0], schemes[1])
    else:
        raise Exception("Path is not of correct length.")


def sub_func(data):
    """
    :param data: list - contains the data needed to run the possible swaps for the bipartition contained in data
    :return: all swapping protocols for the given bipartition.
    """
    scheme_list = []
    G, path, bipartition, symmetric = data
    subpath1, subpath2 = bipartition[0], bipartition[1]
    for schemes in iterate_over_schemes_of_subpaths(subpath1, subpath2, symmetric):
        scheme_list.extend(add_bell_swaps(G, path,
                                          schemes[0], schemes[1]))
    return scheme_list


def get_bipartitions(path, symmetric):
    """from path, get all possible 'bipartitions' that overlap in a single node.
    """
    def path_dict_has_path(path):
        key = return_path_key(path)
        return key in path_dict and path_dict != {}

    path_length = len(path)
    init_length = max(1, math.floor(path_length / 2 - np.log2(path_length)))

    # If uniform_sim, the protocols are symmetric, so swapping between paths of length n1 and n2
    # is equivalent to swapping between paths of length n2 and n1.
    if global_file.params.uniform_sim:
        end_length = math.ceil(path_length/2)
    else:
        end_length = path_length - init_length

    bipartitions = []

    if symmetric:
        subpath1 = path[0:int((path_length+1)/2)]
        subpath2 = path[int((path_length+1)/2-1):]
        bipartitions.append((subpath1, subpath2))
    else:
        for i in range(init_length, end_length):
            subpath1 = path[0:i + 1]
            subpath2 = path[i:path_length]
            if path_dict_has_path(subpath1) and path_dict_has_path(subpath2):
                bipartitions.append((subpath1, subpath2))

    return bipartitions


def iterate_over_schemes_of_subpaths(subpath1, subpath2, symmetric):
    """Returns a generator over all possible combinations of subschemes stored
    in subpath1 and subpath2.
    """

    key1, key2 = return_path_key(subpath1, subpath2)
    if symmetric:
        if key1 == key2:
            return zip(path_dict[key1].values(), path_dict[key2].values())
        else:
            raise ValueError("Keys 1 and 2 for the paths are distinct!")
    else:
        return product(path_dict[key1].values(), path_dict[key2].values())


def iterate_distillation_subschemes(path, general_distillation=False):
    """
    Returns a list of duples of schemes stored over path, which are used for
    distillation protocols. If general_distillation, then we get every possible
    duple of schemes over path. If general_distillation is set to false,
    return only a list of duples of the same scheme. That is,
    [(scheme1, scheme1), (scheme2, scheme2), (scheme3, scheme3), ...]
    A list is returned to prevent the generator from updating as schemes are
    added to path.
    """

    key = return_path_key(path)

    if general_distillation:
        return list(product(path_dict[key].values(), repeat=2))
    else:
        return list(zip(path_dict[key].values(), path_dict[key].values()))


def get_max_length(list_of_paths):
    """ Returns from a list of tuples (paths) the maximum length"""
    return max(map(len, list_of_paths))


def init_elementary_links(G, subpaths, general_distillation):
    """Initialize links (i.e. store all possible schemes) over paths of length 2
    """
    elementary_links, subpaths = get_paths_of_length_k(subpaths, 2)
    if global_file.params.uniform_sim:
        elementary_links = [elementary_links[0]]

    for elementary_link in elementary_links:
        node1, node2 = elementary_link
        # Generate schemes for elementary link generation
        # from both left and right
        if global_file.params.MP_to_IP:
            elementary_link_generation_MP_to_IP(G, node1, node2)
        else:
            if global_file.params.simulation_type == "MP":
                elementary_link_generation_MP(G, node1, node2)
                if not global_file.params.uniform_sim:
                    elementary_link_generation_MP(G, node2, node1)
            elif global_file.params.simulation_type == "IP":
                single_click_protocol(G, node1, node2)
                double_click_protocol(G, node1, node2)
            else:
                print(global_file.params.simulation_type + ' not defined.')

        prune_scheme_from_path(elementary_link)

        distillation_threshold = calc_distillation_threshold(elementary_link)

        for _ in range(global_file.params.distillation_rounds):
            for schemes in iterate_distillation_subschemes(elementary_link,
                                                           general_distillation):
                add_distillation_protocols(
                    G, elementary_link, schemes[0], schemes[1], distillation_threshold)
            prune_scheme_from_path(elementary_link)
    return subpaths
