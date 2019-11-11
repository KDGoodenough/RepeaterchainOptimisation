import numpy as np
import networkx as nx
import random
# from global_file import params
import global_file
import warnings
warnings.simplefilter("ignore", UserWarning)
import matplotlib.pyplot as plt


def create_repeater_chain(distance, n):
    """
    Returns a graph G that is a linear chain of n repeaters
    over the given distance. pert in global_file.py dictates
    whether there should be any random displacements
    """
    G = nx.Graph()
    G.add_node('A', pos=np.array([0, 0]))
    G.add_node('B', pos=np.array([0, distance]))
    G.add_edge('A', 'B')
    if n > 0:
        G = insert_chain(G, 'A', 'B', n)
    nx.set_edge_attributes(G, values=global_file.params.number_of_fibres, name='#fibres')
    return G


def get_ordered_nodes(G, path, subpath1, subpath2):
    """
    Given two subpaths that make up one longer path, returns
    the dictionaries of the nodes at the beginning and end of
    the two subpaths. This is done in the order of begin and
    end nodes of subpath1 and subpath2, depending on their order.
    """
    # Update doc string
    if subpath1 + subpath2 == path:
        node1 = str(subpath1[0])
        node2 = str(subpath1[-1])
        node3 = str(subpath2[0])
        node4 = str(subpath2[-1])
        return G.nodes[node1], G.nodes[node2], G.nodes[node3], G.nodes[node4]

    elif subpath2 + subpath1 == path:
        node1 = str(subpath1[-1])
        node2 = str(subpath1[0])
        node3 = str(subpath2[-1])
        node4 = str(subpath2[0])
        return G.nodes[node1], G.nodes[node2], G.nodes[node3], G.nodes[node4]

    elif subpath1 + subpath2[1:] == path:
        node1 = str(subpath1[0])
        node2 = str(subpath2[0])
        node3 = str(subpath2[-1])
        return G.nodes[node1], G.nodes[node2], G.nodes[node3]

    elif subpath2 + subpath1[1:] == path:
        node1 = str(subpath1[-1])
        node2 = str(subpath1[0])
        node3 = str(subpath2[0])
        return G.nodes[node1], G.nodes[node2], G.nodes[node3]

    else:
        raise ValueError("Subpath1 and subpath2 passed do not form path")


def get_distance_between_nodes(G, path, subpath1, subpath2):
    """
    Returns distance of the link between the two nodes
    connecting subpath1 and subpath2. Returns an error
    if there is no link or if subpath1 and subpath2 do not
    connect to make path.
    """
    if subpath1 + subpath2 == path:
        node1 = str(subpath1[-1])
        node2 = str(subpath2[0])
    elif subpath2 + subpath1 == path:
        node1 = str(subpath2[-1])
        node2 = str(subpath1[0])
    else:
        raise ValueError("Subpath1 and subpath2 passed do not form path")
    return G[node1][node2]['distance']


def get_distance_over_path(G, path):
    """
    Returns distance of the nodes at the beginning and
    ending of path.
    """
    node1 = str(path[0])
    node2 = str(path[-1])

    pos1 = G.nodes[node1]['pos']
    pos2 = G.nodes[node2]['pos']

    return np.sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)


def return_path_key(*args):
    """ Takes a collection of keys (which are tuples of the nodes),
    and returns the key if we have non-uniform simulation, and returns the
    length of that tuple otherwise.
    """
    keys = []
    for arg in args:
        if global_file.params.uniform_sim:
            keys.append(len(arg))
        else:
            keys.append(arg)
    if len(keys) == 1:
        keys = keys[0]
    return keys


def find_candidate_paths(G, Alice, Bob):
    """
    Returns all possible simple paths between Alice and Bob
    (within a cutoff of the shortest path). The paths are returned
    as a list, instead of a generator.

    """
    min_path_length = nx.shortest_path_length(G, source=Alice, target=Bob, weight='weight')
    candidate_paths_gen = nx.all_simple_paths(G, source=Alice, target=Bob)

    candidate_paths = [tuple(candidate_paths)
                       for candidate_paths in candidate_paths_gen]
    return candidate_paths, min_path_length


def find_all_subpaths(all_paths):
    """
    For all paths between the two given end points (given in all_paths),
    find all the possible subpaths. Add them
    to a set to prevent multiple of the same entries
    """
    # Calculate length of the maximum path
    max_length = max(len(s) for s in all_paths)

    subpaths = set()
    for path in all_paths:
        for k in range(0, max_length + 1):
            for ii in range(0, len(path) - k + 1):
                subpaths.add(tuple(path[ii:ii + k]))
    subpaths = filter(None, subpaths)
    return list(subpaths)


def insert_chain(G, A, B, n):
    """
    # This function takes two connected vertices A and B in a graph
    # and replaces the edge with a path graph (a chain of repeaters)
    # of length n. The repeaters are labeled as AB1, AB2,... ABn and
    # are positioned equidistant from each other using the position of
    # A and B.
    """
    # Convert A and B node to str if necessary
    if not isinstance(A, str) or not isinstance(B, str):
        A = str(A)
        B = str(B)

    chain = nx.path_graph(n)  # Create chain of repeaters of length n

    # Get position of A and B, and find the distance between each repeater
    posA = nx.get_node_attributes(G, 'pos')[A]
    posB = nx.get_node_attributes(G, 'pos')[B]
    step = np.subtract(posB, posA)

    # Assign position to each repeater
    repeater_distance = np.round(step / (n + 1), 10)
    try:
        pert = global_file.params.pert
    except AttributeError:
        pert = 0
    for i in chain.nodes():
        posNode = np.add(posA, np.multiply((i + 1), repeater_distance))
        x_offset = random.uniform(-pert, pert)
        y_offset = random.uniform(-pert, pert)
        chain.nodes[i]['pos'] = tuple(np.around(np.add(posNode, (x_offset, y_offset)), 4))
    # Create new chain with nodes AB1, AB2,... ABN
    named_nodes = [A + B + str(i + 1) for i in chain.nodes()]
    # Create mapping from old node names to new names
    mapping = dict(zip(chain.nodes(), named_nodes))
    # Relabel the node names of the chain
    chain = nx.relabel_nodes(chain, mapping)

    # First remove the long edge connecting A and B, and replace by repeater
    # chain
    G.remove_edge(str(A), str(B))
    # Add the edges from the chain including properties to the graph G
    G.add_edges_from(chain.edges(data=True))
    # Add the nodes from the chain including properties to the graph G
    G.add_nodes_from(chain.nodes(data=True))

    # Connect the chain to the graph by connecting A to AB1 and ABN to B
    G.add_edges_from([(str(A), str(A) + str(B) + '1'),
                      (str(B), str(A) + str(B) + str(n))])

    return G


def assign_parameters_to_nodes(G, **kwargs):
    """Assigns arbitrary parameters to all nodes."""
    parameters = get_parameters(**kwargs)
    for key, val in parameters.items():
        nx.set_node_attributes(G, values=val, name=key)
    return G


def assign_parameters_to_specific_nodes(G, nodes, **kwargs):
    """Assigns arbitrary parameters to the list nodes."""
    parameters = get_parameters(**kwargs)
    for node in nodes:
        for key, val in parameters.items():
            G.nodes[node][key] = val
    return G


def get_parameters(**kwargs):
    """
    Returns a subdictionary of the default parameters indicated
    by **kwargs.
    """
    parameters = vars(global_file.params)
    for key, value in kwargs.items():
        parameters[str(key)] = value
    return parameters


def assign_lengths(G):
    """
    Assigns every edge the distance between the connecting vertices
    using the position of the vertices.
    """
    for u, v, d in G.edges(data=True):
        posA = nx.get_node_attributes(G, 'pos')[u]
        posB = nx.get_node_attributes(G, 'pos')[v]

        dist = np.linalg.norm(np.subtract(posA, posB))
        d['distance'] = dist
    return G


def get_paths_of_length_k(subpaths, k):
    """
    From a list of subpaths, returns the subpaths of length k,
    while also removing those subpaths from the list.
    """
    subpaths_of_length_k = [i for i in subpaths if len(
        i) == k]  # all k-length subpaths
    subpaths = [i for i in subpaths if len(i) != k]  # remove k-length subpaths
    return subpaths_of_length_k, subpaths
