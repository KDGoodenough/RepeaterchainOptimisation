# Repeater optimisation

Algorithm for optimising repeater protocols over arbitrary network configurations with arbitrary possible protocols.

## Getting started

### Prerequisites

Requires NetworkX to handle the generation and handling of graphs, graphviz for saving and exporting scheme graphs to .pdf and .dot.

### How it works

#### High-level introduction
The algorithm works almost the same as doi: 10.1073/pnas.0703284104. The input to the function **optimise_schemes**, which performs the actual optimisation, takes as input

- A NetworkX graph G, which represents the network.
- epsilon_p and epsilon_F, which quantifies how fine-grained the optimisation should be. Smaller epsilon results in better results, but longer running times, and vice versa. Negative epsilon will perform no pruning of schemes. See add_scheme(s) defined in scheme_class.py for more detail.
- Alice and Bob, which are the two nodes we want to find the optimal protocol over.

Ignoring pruning for now, the algorithm is as follows:

1. Find all (simple) candidate paths between Alice and Bob over the network and store them in candidate_paths
2. For every candidate_path in candidate_paths, create a list that contains every subpath (without any duplicates).
3. Initialize protocols over `unit length paths'. In other words, for every single node (in subpaths), instantiate parameters of state preparation.
4. For k = 2, for each subpath of length k in subpaths, connect each possible two shorter subpaths of that subpath with every possible connection protocol. This results in new protocols that are stored in path_dict, with the key specified by the tuple of the path. The possible connection protocols are specified and added in add_connection_protocols.
5. For distillation_rounds, for each subpath of length k in subpaths, take two stored protocols, and apply each possible distillation protocol, specified in add_distillation_protocols. Then store each of these protocols in path_dict. Depending on whether general_distillation is set to True or False, take all possible combinations of protocols stored, or only perform distillation over two instances of the same protocol, respectively.
6. remove all subpaths of length k from subpaths, and increase set k to k + 1. Repeat the above until subpaths is empty.


The amount of possible protocols can grow very quickly. In the paper we discuss the various ways we make sure the number of considered schemes remains tractable.

The meat of the algorithm is in optimise_schemes, which takes as inputs:

    :param G: nx.Graph. Graph over which to perform the optimisation. The nodes of G should have
        several attributes, pos for absolute position of the node, and error parameters.
    :param general_distillation: bool - whether to perform distillation over two of the
        same pair (=False) or any possible combination (=True)
    :param Alice: ID used for Alice
    :param Bob: ID used for Bob
    :param symmetricoptimisation: bool - indicates whether the bisection heuristic is to be applied or not
    """

The parameters can be uniformly changed by changing them in global_file.py. Changing parameters non-uniformly can be done by using the functions in network_functions.

The optimisation results are contained in the dictionary path_dict, where the keys are given by the paths, and the values by the class Schemes.
Schemes is a dictionary, with the keys given by rounded fidelity and probability tuples, and the values by schemes.
Thus, the optimisation results for Alice and Bob can be found in path_dict[path_between_alice_and_bob].

To add your own protocols, add them to either connection_protocols or distillation_protocols:

    For ELG, add them to init_elementary_links in connection_and_distillation
    For swapping, add them to add_bell_swaps in connection_protocols
    for distillation, add them to add_distillation_protocols in distillation_protocols

A protocol at the very least needs to call at the end add_scheme with a scheme. The scheme contains:

    :param path: hashable. path over which the entanglement has been generated
    :param protocoltype: str. label for the last protocol used
    :param subscheme1: Scheme or None. One of the schemes used to generate the scheme. None if scheme is associated to ELG
    :param subscheme2: Scheme or None. One of the schemes used to generate the scheme. None if scheme is associated to ELG
    :param time: float. Float indicating the generation time for the scheme
    :param state: matrix. Matrix representing the average state generated with the scheme
    :param prob: float. probability that the scheme has succeeded
    :param fidelity: float. fidelity of the state to the maximally entangled state
    
An example code can be found in main.py, where the optimisation is run for information-processig platforms,
for a distance of 20 km, n=1 repeater nodes, and the parameters in global_file.py. The results are displayed in a
fidelity vs time plot.
# RepeaterOptimisationCode
