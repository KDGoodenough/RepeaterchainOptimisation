from global_file import path_dict
import collections
from subprocess import check_call
import networkx as nx
import numpy as np
from numpy.ma import masked_invalid
import qutip as qt
from networkx.drawing.nx_pydot import write_dot
from network_functions import get_distance_over_path
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import pickle
import os
from pathlib import Path
import glob
from scipy.stats import entropy
import scipy.ndimage as ndimage
from random import randint
import matplotlib
import re
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as ticker
rcParams.update({'figure.autolayout': True})
matplotlib.rc('text', usetex=True)


def collect_optimal_protocols(G, Alice, Bob, stored_path_dict=None, n=None, uniform_sim=True):
    # Path_dict is passed in case the data was stored on file
    # Might be moved to scheme_class
    """
    From all candidate_paths (which should be paths
    between Alice and Bob) return (ordered) optimal
    fidelities, times and schemes.
    """

    # Retrieves all schemes stored between Alice and B, unless an n
    # is specified, which is the length over which to collect optimal protocols while
    # running uniform simulations.
    if uniform_sim:
        if G:
            scheme_list = path_dict[G.number_of_nodes()].values()
        else:
            if not n:
                n = max(stored_path_dict)
            scheme_list = stored_path_dict[n].values()
    else:
        candidate_paths = list(nx.all_simple_paths(G, source=Alice, target=Bob))
        scheme_list = [path_dict[tuple(path)].values() for path in candidate_paths]
        scheme_list = [
            scheme for schemes in scheme_list for scheme in schemes]

    ordered_schemes = sorted(scheme_list, key=lambda x: x.fidelity, reverse=True)
    optimal_schemes = []
    for index, scheme in enumerate(ordered_schemes):
        if index == 0:
            optimal_schemes.append(scheme)
            previous_time = scheme.time
        else:
            if scheme.fidelity == optimal_schemes[0].fidelity:
                if scheme.time < previous_time:
                    optimal_schemes[0] = scheme
                    previous_time = scheme.time
            else:
                if scheme.time < previous_time:
                    optimal_schemes.append(scheme)
                    previous_time = scheme.time

    optimal_schemes.reverse()

    fidelities = [scheme.fidelity for scheme in optimal_schemes]
    generation_times = [scheme.time for scheme in optimal_schemes]
    return fidelities, generation_times, optimal_schemes


def plot_max_fidelity_with_time(G=None, Alice=None, Bob=None, stored_path_dict=None, n=None, color=None):
    fidelities, opt_times, _ = collect_optimal_protocols(G, Alice, Bob, stored_path_dict, n)
    plt.scatter(fidelities[-1], opt_times[-1], c=color, s=57)
    L = get_distance_over_path(G, (Alice, Bob))
    if int(L) != 15:
        plt.annotate('L = ' + str(int(L)), (fidelities[-1]+0.0025, opt_times[-1]+0.00001), fontsize=16)
    else:
        plt.annotate('L = ' + str(int(L)), (fidelities[-1]+0, opt_times[-1]+0.00005), fontsize=16)


def plot_fidelity_vs_time_plot(G=None, Alice=None, Bob=None, stored_path_dict=None, threshold_fidelity=0.5,
                               n=None, label='No label passed', color=None, connected=False, linestyle='-', ax=None):
    """
    For every scheme between Alice and Bob, generate a plot
    with the optimal fidelities and times. label is the label for the legend
    connected is used to indicate whether the graph should be connected.
    """

    fidelities, opt_times, _ = collect_optimal_protocols(G, Alice, Bob, stored_path_dict, n)

    # Tries to prepend trivially easy protocols.
    # Might not work if opt_times is empty
    try:
        opt_times.insert(0, opt_times[0])
        fidelities.insert(0, threshold_fidelity)
    except BaseException:
        pass

    lw = 2
    alpha = 1
    if connected:
        if color:
            if ax:
                ax.step(fidelities, opt_times, linestyle=linestyle, where='pre',
                        linewidth=lw, label=label, ms=18, color=color, alpha=alpha)
            else:
                plt.step(fidelities, opt_times, linestyle=linestyle, where='pre',
                         linewidth=lw, label=label, ms=18, color=color, alpha=alpha)
        else:
            if ax:
                ax.step(fidelities, opt_times, linestyle=linestyle, where='pre',
                        linewidth=lw, label=label, ms=18, alpha=alpha, ax=ax)
            else:
                plt.step(fidelities, opt_times, linestyle=linestyle, where='pre',
                         linewidth=lw, label=label, ms=18, alpha=alpha)
    else:
        plt.plot(fidelities, opt_times, 'o', label=label, ms=18)

    log_plot = True
    if log_plot:
        if ax:
            ax.set_yscale('log')
        else:
            plt.yscale('log')

    if ax:
        ax.tick_params(axis='both', labelsize=14)
        ax.set_xticks(np.arange(0.5, 1.1, step=0.1))
    else:
        plt.xticks(fontsize=28)
        plt.yticks(fontsize=28)
    plt.xlim(threshold_fidelity, 1)


def generate_pdf_of_scheme_graph(scheme, node_ID, general_distillation=False):
    """
    For scheme, generate, save and open the corresponding scheme graph with given filename.
    general_distillation set to False implies we only do repeated distillation, saving
    reducing clutter in the scheme graph.
    """
    scheme_graph = nx.DiGraph()
    scheme_graph = create_scheme_tree(scheme, scheme_graph, general_distillation, node_ID=node_ID)
    return scheme_graph


def create_scheme_tree(schemes, scheme_graph, general_distillation, node_ID):
    for scheme in get_iterable(schemes):
        if scheme:
            node_ID1 = randint(0, 999999)
            node_ID2 = randint(0, 999999)
            if scheme.subscheme1.protocoltype[0:2] not in ['Ns', 'sc', 'dc']:
                scheme_graph = add_protocol_to_graph_symmetric(scheme_graph, scheme,
                                                               scheme.subpath1, general_distillation,
                                                               node_ID, node_ID1)
                create_scheme_tree(scheme.subscheme1,
                                   scheme_graph, general_distillation, node_ID1)
            else:
                label_parent = simplify_labels(scheme.protocoltype)
                label_child = simplify_labels(scheme.subscheme1.protocoltype)
                scheme_graph.add_node(node_for_adding=node_ID, node_size=300, node_shape='o')
                scheme_graph.add_node(node_for_adding=node_ID1, label='2', node_size=300, node_shape='o')
                node_ID3 = randint(0, 99999999)
                scheme_graph.add_node(node_for_adding=node_ID3, label='ELG', node_size=400, node_shape='o')
                scheme_graph.add_edge(node_ID, node_ID1, label=label_parent)
                scheme_graph.add_edge(node_ID1, node_ID3, label=label_child)

            if scheme.subscheme2.protocoltype[0:2] not in ['Ns', 'sc', 'dc']:
                scheme_graph = add_protocol_to_graph_symmetric(scheme_graph, scheme,
                                                               scheme.subpath2, general_distillation,
                                                               node_ID, node_ID2)
                create_scheme_tree(scheme.subscheme2,
                                   scheme_graph, general_distillation, node_ID2)
            else:
                label_parent = simplify_labels(scheme.protocoltype)
                label_child = simplify_labels(scheme.subscheme2.protocoltype)

                scheme_graph.add_node(node_for_adding=node_ID, node_size=300, node_shape='o')
                scheme_graph.add_node(node_for_adding=node_ID2, label='2', node_size=300, node_shape='o')
                node_ID3 = randint(0, 99999999)
                scheme_graph.add_node(node_for_adding=node_ID3, label='ELG', node_size=400, node_shape='o')
                scheme_graph.add_edge(node_ID, node_ID2, label=label_parent)
                scheme_graph.add_edge(node_ID2, node_ID3, label=label_child)

    return scheme_graph


def add_protocol_to_graph(scheme_graph, scheme, subpath, general_distillation):
    """ Add the protocol stored in scheme to scheme_graph. """
    # Get simple path names that 'connect' the protocols and add them to
    # scheme.graph.
    path_name = get_simple_path_name(scheme.path)
    subpath_name = get_simple_path_name(subpath)
    scheme_graph.add_node(path_name)
    scheme_graph.add_node(subpath_name)
    # If general_distillation = False, we can simplify the notation used in
    # the scheme_graph significantly
    if general_distillation:
        if scheme_graph.has_edge(path_name, subpath_name):
            label = scheme_graph.get_edge_data(path_name, subpath_name)[
                'label'] + "-" + scheme.protocoltype
        else:
            label = scheme.protocoltype
        scheme_graph.add_edge(path_name, subpath_name, label=label)
    else:
        if scheme_graph.has_edge(path_name, subpath_name):
            iterations = scheme_graph.get_edge_data(
                path_name, subpath_name)['iterations'] + 1
            label = str(iterations) + ' x ' + scheme.protocoltype
            scheme_graph.add_edge(path_name, subpath_name,
                                  label=label, iterations=iterations)
        else:
            label = scheme.protocoltype
            scheme_graph.add_edge(path_name, subpath_name,
                                  label=label, iterations=1)
    return scheme_graph


def add_protocol_to_graph_symmetric(scheme_graph, scheme, subpath, general_distillation, node_ID_parent, node_ID_child):
    """ Add the protocol stored in scheme to scheme_graph. """
    # Get simple path names that 'connect' the protocols and add them to scheme.graph.

    path_name = str(len(scheme.path))
    if subpath:
        subpath_name = str(len(subpath))
    else:
        subpath_name = str(1)

    scheme_graph.add_node(node_for_adding=node_ID_parent, label=path_name, node_size=300, node_shape='o')
    if subpath:
        scheme_graph.add_node(node_for_adding=node_ID_child, label=subpath_name, node_size=300, node_shape='o')

        if general_distillation:
            if scheme_graph.has_edge(path_name, subpath_name):
                raise ValueError("There cannot be an edge between path_name"
                                 " and subpath_name since the vertices are new.")
            else:
                label = simplify_labels(scheme.protocoltype)
            scheme_graph.add_edge(node_ID_parent, node_ID_child, label=label)
        else:
            raise ValueError("General distillation needs to be true for symmetric scheme representation.")
        return scheme_graph
    else:
        label = simplify_labels(scheme.protocoltype)
        scheme_graph.add_node(node_for_adding=node_ID_child, label=label, node_size=300, node_shape='o')
        scheme_graph.add_edge(node_ID_parent, node_ID_child)
        return scheme_graph


def simplify_labels(label):
    if label[0:2] == "Ns":
        floats = re.findall(r"[-+]?\d*\.\d+|\d+", label)
        a = float(floats[0])
        b = int(floats[1])
        return "Ns = {:3.4f},\nr = {}".format(a, b)
    elif label[0:5] == "cSwap":
        floats = re.findall(r"[-+]?\d*\.\d+|\d+", label)
        assert(len(floats) == 1)
        return "Swap,\nr = " + str(floats[0])
    elif label[0:2] == "dc":
        floats = re.findall(r"[-+]?\d*\.\d+|\d+", label)
        assert (len(floats) == 1)
        return "DC,\nr = " + str(floats[0])
    elif label[0] == "d":
        floats = re.findall(r"[-+]?\d*\.\d+|\d+", label)
        assert (len(floats) == 1)
        return "Distill,\nr = " + str(floats[0])
    elif label[0:2] == "sc":
        floats = re.findall(r"[-+]?\d*\.\d+|\d+", label)
        a = float(floats[0])
        b = int(floats[1])
        return "$\\theta$ = {:3.4f},\nr = {}".format(a, b)
    else:
        return label


def get_simple_path_name(path_name):
    """ Simplifies path name by removing brackets, hyphens and spaces """
    simple_path_name = str(path_name).replace("(", "")
    simple_path_name = simple_path_name.replace(")", "")
    simple_path_name = simple_path_name.replace("'", "")
    simple_path_name = simple_path_name.replace(",", "-")
    simple_path_name = simple_path_name.replace(" ", "")
    if simple_path_name[-1] == "-":
        simple_path_name = simple_path_name[0:-1]
    return simple_path_name


def plot_scheme_density(path):
    schemes = path_dict[path]
    fidelities_rounded = []
    probabilities_rounded = []

    fidelities = []
    probabilities = []
    times = []

    for key in schemes:
        fidelities_rounded.append(key[0])
        probabilities_rounded.append(key[1])

        fidelities.append(schemes[key].fidelity)
        probabilities.append(schemes[key].prob)
        times.append(schemes[key].time)

    cm = plt.cm.get_cmap('OrRd')

    norm_times = [(time - min(times)) / (max(times) - min(times))
                  for time in times]
    plt.scatter(
        fidelities_rounded,
        probabilities_rounded,
        c=norm_times,
        cmap=cm)
    plt.colorbar()
    plt.figure()
    plt.scatter(fidelities, probabilities, c=norm_times, cmap=cm)
    plt.colorbar()
    plt.show()


def select_and_generate_plots(colors=None, labels=None, linestyles=None, ax=None):
    files = []
    filename = "init_filename"
    while filename:
        Tk().withdraw()  # we don't want a full GUI, so keep the root window from appearing
        filename = askopenfilename()  # show an "Open" dialog box and return the path to the selected file
        files.append(filename)

    if not colors:
        colors = ['r', 'g', 'b', 'y', 'c', 'k', 'r', 'g', 'b', 'y', 'c']
    else:
        colors = colors*10

    if labels is not None:
        labels = labels*10

    for index, file in enumerate(files):
        if file:

            if linestyles is not None:
                linestyle = linestyles[index]
            else:
                linestyle = '-'

            if not labels:
                label = Path(file).stem
            else:
                label = labels[index]
            with open(file, 'rb') as f:
                try:
                    stored_path_dict = pickle.load(f)
                    if stored_path_dict:
                        plot_fidelity_vs_time_plot(None, None, None, stored_path_dict=stored_path_dict,
                                                   label=label, color=colors[index],
                                                   connected=True, linestyle=linestyle, ax=ax)
                except EOFError:
                    pass


def select_and_generate_fixed_interrepeater_distance_plots(colors=None, labels=None, linestyles=None, ax=None):
    Tk().withdraw()  # we don't want a full GUI, so keep the root window from appearing
    file = askopenfilename()  # show an "Open" dialog box and return the path to the selected file

    if not colors:
        colors = ['r', 'g', 'b', 'y', 'c', 'k', 'r', 'g', 'b', 'y', 'c']
    else:
        colors = colors*10

    if labels is not None:
        labels = labels*10

    if file:
        with open(file, 'rb') as f:
            stored_path_dict = pickle.load(f)
            for i in range(2, max(stored_path_dict) + 1):
                if not labels:
                    label = Path(file).stem
                else:
                    label = labels[i-2]
                if linestyles is not None:
                    linestyle = linestyles[i-2]
                else:
                    linestyle = '-'

                plot_fidelity_vs_time_plot(None, None, None, stored_path_dict=stored_path_dict, n=i,
                                           label=Path(file).stem, color=colors[i-2], connected=True,
                                           linestyle=linestyle, ax=ax)


def select_data():
    Tk().withdraw()
    file = askopenfilename()
    if file:
        with open(file, 'rb') as f:
            return pickle.load(f)


def plot_connection_density(path):
    schemes = path_dict[path]
    for scheme in schemes:
        if scheme.protocoltype[0] == 'c':
            f1, f2 = scheme.subscheme1.fidelity, scheme.subscheme2.fidelity
            p1, p2 = scheme.subscheme1.prob, scheme.subscheme2.prob
            plt.plot([f1, f2], [p1, p2])

    plt.show()


def save_and_open_pdf_of_scheme_graph(graph, filename):
    """ Writes tree graph to pdf and dot file with given filename """
    write_dot(graph, filename + ".dot")
    os.environ['PATH'] += ":" + "/usr/local/bin"
    check_call(["dot", "-T", "pdf", filename +
                ".dot", "-o", filename + ".pdf"])
    os.system("open " + filename + ".pdf")


def get_iterable(x):
    """
    This function can be used when we want to iterate over the list x,R
    but x might also be a singleton.
    """
    if isinstance(x, collections.Iterable):
        return x
    else:
        return (x,)


def get_name_of_path(n, endpoint=False):
    # This is only used for testing purposes, and returns the
    # tuple ('A', 'AB1', 'AB2', 'AB3',...'ABn', 'B')

    path_name = ['A']
    for i in range(0, n):
        path_name += ['AB' + str(i+1)]
    if endpoint:
        path_name += ['B']
    return tuple(path_name)


def gather_schemes_from_files():
    files = []
    filename = "init_filename"
    while filename:
        Tk().withdraw()  # we don't want a full GUI, so keep the root window from appearing
        filename = askopenfilename()  # show an "Open" dialog box and return the path to the selected file
        files.append(filename)

    fidelities_list = []
    times_list = []
    for file in files:
        if file:
            with open(file, 'rb') as f:
                stored_path_dict = pickle.load(f)
                fidelities, generation_times, _ = collect_optimal_protocols(None, "A", "B",
                                                                            stored_path_dict=stored_path_dict, n=None)
                fidelities_list.append(fidelities)
                times_list.append(generation_times)

    return fidelities_list, times_list


def generate_param_colour_plot(F, ax, scale):
    def fmt(x):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    text = "15 km, 0 rptrs * NP, T=1e-2-1e1, P = 0.9-0.95, march 6 *.pkl"

    x = 20
    y = 20

    files = glob.glob(text)
    files.sort(key=lambda x: [int(s) for s in x.split() if s.isdigit()][2])

    rates = [[0 for i in range(x)] for j in range(y)]

    counter1 = 0
    for file in files:
        with open(file, 'rb') as f:
            stored_path_dict = pickle.load(f)
            if stored_path_dict:
                fidelities, opt_times, opt_schemes = collect_optimal_protocols(G=None, Alice=None, Bob=None,
                                                                               stored_path_dict=stored_path_dict)
                try:
                    # index of the smallest fidelity greater than F
                    index = next(x[0] for x in enumerate(fidelities) if x[1] >= F)
                    max_rate = 1/opt_times[index]
                except StopIteration:  # exception occurs when there is no scheme that achieves F, return infinite time
                    max_rate = 0
                counter = [int(s) for s in file.split() if s.isdigit()][2] - 1
                print(counter)
                print(file)

                rates[counter % y][counter // y] = max_rate

        counter1 += 1

    plt.set_cmap('RdPu')
    current_cmap = matplotlib.cm.get_cmap()
    current_cmap.set_bad(color='black')
    try:
        min_time = np.nanmin(np.nanmin(rates))
        max_time = np.nanmax(np.nanmax(rates))
    except:
        min_time, max_time = 0, 0
    times = masked_invalid(rates)

    current_cmap = matplotlib.cm.get_cmap()
    current_cmap.set_bad(color='black', alpha=1.)

    skr_grad = np.gradient(ndimage.gaussian_filter(times, sigma=0.35, order=0))
    lengths = [[np.power(0.00001+skr_grad[0][j][i]**2 + skr_grad[1][j][i]**2, 1/2)
                if np.isfinite(skr_grad[0][j][i]) or np.isfinite(skr_grad[1][j][i])
                else 1 for i in range(x)] for j in range(y)]

    skr_grad[0] = skr_grad[0]/lengths
    skr_grad[1] = skr_grad[1]/lengths

    ax.quiver(skr_grad[1], -skr_grad[0], scale=scale, alpha=0.8)

    im = ax.matshow(rates, cmap=current_cmap, interpolation=None)

    t = np.linspace(min_time, max_time, 5)
    cb = plt.colorbar(im, ax=ax, ticks=t, format=ticker.FuncFormatter(fmt), fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=11)

    cb.ax.set_title('Generation rate', fontsize=10)


def generate_param_colour_plot_throughput(ax):
    def fmt(x):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    def binary_entropy(x):
        if x not in [0, 1]:
            return -x*np.log2(x)-(1-x)*np.log2(1-x)
        else:
            return 0

    def calc_throughput(scheme):
        state = scheme.state
        time = scheme.time
        simulation_type = "IP"
        if simulation_type == "IP":
            proj1 = qt.Qobj([[1, 0, 0, 0],
                             [0, 0, 0, 0],
                             [0, 0, 0, 0],
                             [0, 0, 0, 1]])
            proj2 = (1/2)*qt.Qobj([[1, 0, 0, 1],
                                   [0, 1, 1, 0],
                                   [0, 1, 1, 0],
                                   [1, 0, 0, 1]])
            proj3 = (1/2)*qt.Qobj([[1, 0, 0, -1],
                                   [0, 1, 1, 0],
                                   [0, 1, 1, 0],
                                   [-1, 0, 0, 1]])

            state00 = qt.Qobj([[1/2, 0, 0, 1/2],
                               [0,   0, 0, 0],
                               [0,   0, 0, 0],
                               [1/2, 0, 0, 1/2]])

            state01 = qt.Qobj([[1/2, 0, 0, -1/2],
                               [0,   0, 0, 0],
                               [0,   0, 0, 0],
                               [-1/2, 0, 0, 1/2]])

            state10 = (1/2)*qt.Qobj([[0, 0, 0, 0],
                                     [0, 1, 1, 0],
                                     [0, 1, 1, 0],
                                     [0, 0, 0, 0]])

            state11 = (1/2)*qt.Qobj([[0, 0, 0, 0],
                                     [0, 1, -1, 0],
                                     [0, -1, 1, 0],
                                     [0, 0, 0, 0]])

            proj1.dims = [[2, 2], [2, 2]]
            proj2.dims = [[2, 2], [2, 2]]
            proj3.dims = [[2, 2], [2, 2]]
            state00.dims = [[2, 2], [2, 2]]
            state01.dims = [[2, 2], [2, 2]]
            state10.dims = [[2, 2], [2, 2]]
            state11.dims = [[2, 2], [2, 2]]

        elif simulation_type == "MP":
            proj1 = qt.ket2dm(qt.basis(81, 28)) + qt.ket2dm(qt.basis(81, 12))
            proj2 = qt.ket2dm(qt.basis(81, 28) + qt.basis(81, 12)) / 2 + \
                    qt.ket2dm(qt.basis(81, 28) - qt.basis(81, 12)) / 2
            proj3 = qt.ket2dm(qt.basis(81, 28) + 1j*qt.basis(81, 12)) / 2 + \
                    qt.ket2dm(qt.basis(81, 28) - 1j*qt.basis(81, 12)) / 2

            vec01 = qt.basis(9, 1)
            vec10 = qt.basis(9, 3)

            state00 = qt.ket2dm(qt.tensor(vec01, vec10) + qt.tensor(vec10, vec01))/2
            state01 = qt.ket2dm(qt.tensor(vec01, vec10) - qt.tensor(vec10, vec01))/2
            state10 = qt.ket2dm(qt.tensor(vec01, vec01) + qt.tensor(vec10, vec10))/2
            state11 = qt.ket2dm(qt.tensor(vec01, vec01) - qt.tensor(vec10, vec10))/2

            proj1.dims = [[9, 9], [9, 9]]
            proj2.dims = [[9, 9], [9, 9]]
            proj3.dims = [[9, 9], [9, 9]]
            state00.dims = [[9, 9], [9, 9]]
            state01.dims = [[9, 9], [9, 9]]
            state10.dims = [[9, 9], [9, 9]]
            state11.dims = [[9, 9], [9, 9]]

        bb84 = False
        if bb84:
            qber1 = (state * proj1).tr()
            qber2 = (state * proj2).tr()
            skf = max(1 - binary_entropy(qber1) - binary_entropy(qber2), 0)  # calc secret_key fraction
        else:
            p00 = (state * state00).tr()
            p01 = (state * state01).tr()
            p10 = (state * state10).tr()
            p11 = (state * state11).tr()

            normalisation = p00+p01+p10+p11

            # define P' distribution as in https://iopscience.iop.org/article/10.1088/2058-9565/aab31b, eq D1 to D8
            Px0 = (p00+p01)**2 + (p10+p11)**2
            Px1 = 1-Px0
            p00n = (p00**2+p01**2)/Px0
            p01n = (2*p00*p01)/Px0
            p10n = (p10**2+p11**2)/Px0
            p11n = (2*p10*p11)/Px0
            if ((p00 + p01) * (p10 + p11)) == 0:
                Pbinary_entropy = 1
            else:
                Pbinary_entropy = (p00 * p10 + p01 * p11) / ((p00 + p01) * (p10 + p11))

            skf = max(1-entropy([p00, p01, p10, p11], base=2)+(Px1/2)*binary_entropy(Pbinary_entropy),
                      (Px0/2)*(1-entropy([p00n, p01n, p10n, p11n], base=2)), 0)

        return skf/time

    text = "50 km, 1 rptrs * G=0.98-1 in 20, P=* in 20, eps = 0.01, IP, 1 node, 50km, 8 feb, *.pkl"
    x = 20
    y = 20

    files = glob.glob(text)
    files = [file for file in files if "pessimistic" not in file]
    files.sort(key=lambda x: [int(s) for s in x.split() if s.isdigit()][2])
    files.sort(key=lambda x: [int(s) for s in x.split() if s.isdigit()][1])
    skr = [[0 for i in range(x)] for j in range(y)]
    for file in files:
        print(file)
    counter = 0
    for file in files:
        with open(file, 'rb') as f:
            stored_path_dict = pickle.load(f)
            if stored_path_dict:
                fidelities, opt_times, opt_schemes = collect_optimal_protocols(G=None, Alice=None, Bob=None,
                                                                               stored_path_dict=stored_path_dict)
                throughputs = [calc_throughput(scheme) for scheme in opt_schemes]
                max_throughput = max(throughputs)
                counter = [int(s) for s in file.split() if s.isdigit()][2] - 1
                skr[counter % y][counter // y] = max_throughput

        counter += 1

    print(counter)
    plt.set_cmap('RdPu')
    current_cmap = matplotlib.cm.get_cmap()
    current_cmap.set_bad(color='black')
    try:
        min_time = np.nanmin(np.nanmin(skr))
        max_time = np.nanmax(np.nanmax(skr))
    except:
        min_time, max_time = 0, 0

    skr = masked_invalid(skr)

    current_cmap = matplotlib.cm.get_cmap()
    current_cmap.set_bad(color='black', alpha=1.)

    skr_grad = np.gradient(skr)
    lengths = [[np.power(0.001+skr_grad[0][j][i]**2 + skr_grad[1][j][i]**2, 1/2)
                if np.isfinite(skr_grad[0][j][i]) or np.isfinite(skr_grad[1][j][i])
                else 1 for i in range(x)] for j in range(y)]

    maxlength = max(max(lengths))*1.5
    skr_grad[0] = skr_grad[0]/maxlength
    skr_grad[1] = skr_grad[1]/maxlength
    plt.quiver(skr_grad[1], -skr_grad[0])

    im = ax.matshow(skr, cmap=current_cmap, interpolation=None)

    plt.rc('text', usetex=True)
    t = np.linspace(min_time, max_time, 5)
    cb = plt.colorbar(im, ax=ax, ticks=t, format=ticker.FuncFormatter(fmt), fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=16)
    cb.ax.set_title('Secret-key rate', fontsize=17, pad=12)
