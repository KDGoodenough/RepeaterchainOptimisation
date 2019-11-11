
from ete3 import Tree, TreeNode
from plot_functions import generate_pdf_of_scheme_graph, select_data
from plot_functions import save_and_open_pdf_of_scheme_graph
import networkx as nx
from random import randint
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import matplotlib
from matplotlib.figure import figaspect
matplotlib.use('TkAgg')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt
import nxpd
from math import atan2, degrees


def get_label(G, n):
    k = G.predecessors(n)
    labels = []
    for i in k:
        labels.append(G[i][n]['label'])
    if len(labels) == 2:
        assert labels[0] == labels[1]
        return labels[0]
    elif len(labels) == 0:
        return ""
    elif len(labels) == 1:
        return labels[0]
    return labels[0]


def get_label_pos(G, n, pos):
    k = G.predecessors(n)
    pos2 = []
    for i in k:
        pos2.append(pos[i])

    if len(pos2) == 2:
        return ((pos2[0][0] + pos2[1][0] + pos[n][0])/3, (pos2[0][1] + pos2[1][1] + pos[n][1])/3)
    elif len(pos2) == 0:
        return (pos[n][0], pos[n][1]-10)
    elif len(pos2) == 1:
        return ((pos2[0][0] + pos[n][0])/2, (pos2[0][1] + pos[n][1])/2)

# G = nx.DiGraph()
#
#
#
# G.add_node(1, node_size=1800, label='$\\textrm{QR}_{i}$-$\\textrm{QR}_{j}$')
# G.add_node(2, node_size=1800, label='$\\textrm{QR}_{j}$-$\\textrm{QR}_{k}$')
# G.add_node(3, node_size=1800, label='$\\textrm{QR}_{i}$-$\\textrm{QR}_{k}$')
# G.add_edge(3,1, label='Swap,\nr = r$^*$')
# G.add_edge(3,2, label='Swap,\nr = r$^*$')
# #save_and_open_pdf_of_scheme_graph(G, 'test11june.png')
# #pos = graphviz_layout(G, prog='dot')
# labels = nx.get_node_attributes(G, 'label')
# #print(pos)
# edge_labels = nx.get_edge_attributes(G, 'label')
# node_sizes = [size for (_, size) in nx.get_node_attributes(G, 'node_size').items()]
# # node_shapes = [shape for (_, shape) in nx.get_node_attributes(G, 'node_shape').items()]
# #node_shapes = nx.get_node_attributes(G, 'node_shape')
#
# G = G.reverse()
# #f = plt.figure()
# # width = 6
# width = 3
# height = 2
# plt.figure(figsize=(width, height))
# pos = {1: (240, 25.0), 2: (260, 25.0), 3: (250, 54.0)}
# # pos = {1: (5, 5), 2: (5, 5.0), 3: (5, 5)}
# nx.draw(G, pos, with_labels=False, arrows=True, node_color='orchid', node_size=node_sizes,font_size=4)
# # nx.draw_networkx_labels(G, pos, with_labels=False, arrows=True, node_color='orchid',
# #                         node_size=node_sizes,font_size=8)
#
# text = nx.draw_networkx_labels(G, pos=pos, font_size=8, labels=labels, rotation='horizontal')
# # for i, t in text.items():
# #     t.set_rotation(90)
#
# # nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, bbox=dict(alpha=0), font_size=7)
# print(pos)
# label_pos = {node: get_label_pos(G, node, pos) for (node, position) in pos.items()}
# labels = {node: get_label(G, node) for (node, position) in pos.items()}
# # print(label_pos)
# # label_pos = {1: get_label_pos(G, 1, pos)}
# text = nx.draw_networkx_labels(G, pos=label_pos, font_size=12, labels=labels, rotation='horizontal')
# for i, t in text.items():
#     if t.get_text()[0:2] in ['Ns', 'sc', 'dc', "$\\", r"$\""]:
#         for k in G.predecessors(i):
#             diff = (label_pos[i][0]-label_pos[k][0], label_pos[i][1]-label_pos[k][1])
#             angle = degrees(atan2(diff[1], diff[0]))
#             print(angle)
#         t.set_rotation('vertical')
#     else:
#         t.set_rotation('vertical')
#         # pass
# #nxpd.draw(G)
# # plt.show()
# #
# fig = plt.gcf()
# ax = plt.gca()
# print(ax.get_xlim())
# print(ax.get_ylim())
# ax.set_xlim((235, 265))
# ax.set_ylim((18, 64))
# size = fig.get_size_inches()
# # plt.gcf().su#bplots_adjust(bottom=0.15)
# plt.savefig('Notes/scheme_graphs/swap_graph.pdf')#, bbox_inches='tight')
#
#
#
#
#
# G = nx.DiGraph()
#
#
#
# G.add_node(1, node_size=2100, label='ELG')
# G.add_node(2, node_size=2000, label='$\\textrm{QR}_{i}$-$\\textrm{QR}_{i+1}$')
# G.add_edge(2,1, label='ELG,\nr = r$^*$')
# #save_and_open_pdf_of_scheme_graph(G, 'test11june.png')
# pos = graphviz_layout(G, prog='dot')
# labels = nx.get_node_attributes(G, 'label')
# #print(pos)
# edge_labels = nx.get_edge_attributes(G, 'label')
# node_sizes = [size for (_, size) in nx.get_node_attributes(G, 'node_size').items()]
# # node_shapes = [shape for (_, shape) in nx.get_node_attributes(G, 'node_shape').items()]
# #node_shapes = nx.get_node_attributes(G, 'node_shape')
#
# G = G.reverse()
# #f = plt.figure()
# # width = 6
# width = 3
# height = 2
# plt.figure(figsize=(width, height))
# #pos = {1: (240, 25.0), 2: (260, 25.0), 3: (250, 54.0)}
# # pos = {1: (5, 5), 2: (5, 5.0), 3: (5, 5)}
# nx.draw(G, pos, with_labels=False, arrows=True, node_color='orchid', node_size=node_sizes,font_size=4)
# # nx.draw_networkx_labels(G, pos, with_labels=False, arrows=True, node_color='orchid',
# #                         node_size=node_sizes,font_size=8)
#
# text = nx.draw_networkx_labels(G, pos=pos, font_size=8, labels=labels, rotation='horizontal')
# # for i, t in text.items():
# #     t.set_rotation(90)
#
# # nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, bbox=dict(alpha=0), font_size=7)
# print(pos)
# label_pos = {node: get_label_pos(G, node, pos) for (node, position) in pos.items()}
# labels = {node: get_label(G, node) for (node, position) in pos.items()}
# # print(label_pos)
# # label_pos = {1: get_label_pos(G, 1, pos)}
# text = nx.draw_networkx_labels(G, pos=label_pos, font_size=12, labels=labels, rotation='horizontal')
# for i, t in text.items():
#     if t.get_text()[0:2] in ['Ns', 'sc', 'dc', "$\\", r"$\""]:
#         for k in G.predecessors(i):
#             diff = (label_pos[i][0]-label_pos[k][0], label_pos[i][1]-label_pos[k][1])
#             angle = degrees(atan2(diff[1], diff[0]))
#             print(angle)
#         t.set_rotation('vertical')
#     else:
#         t.set_rotation('vertical')
#         # pass
# #nxpd.draw(G)
# # plt.show()
#
# fig = plt.gcf()
# ax = plt.gca()
# print(ax.get_xlim())
# print(ax.get_ylim())
# # ax.set_xlim((235, 265))
# ax.set_ylim((-8, 142))
# size = fig.get_size_inches()
# # plt.gcf().su#bplots_adjust(bottom=0.15)
# plt.savefig('Notes/scheme_graphs/ELG_graph.pdf')#, bbox_inches='tight')
#
#
#
#
#
#
# G = nx.DiGraph()
#
#
#
# G.add_node(1, node_size=1800, label='$\\textrm{QR}_{i}$-$\\textrm{QR}_{j}$')
# G.add_node(2, node_size=1800, label='$\\textrm{QR}_{i}$-$\\textrm{QR}_{j}$')
# G.add_node(3, node_size=1800, label='$\\textrm{QR}_{i}$-$\\textrm{QR}_{j}$')
# G.add_edge(3,1, label='Distill,\nr = r$^*$')
# G.add_edge(3,2, label='Distill,\nr = r$^*$')
# #save_and_open_pdf_of_scheme_graph(G, 'test11june.png')
# #pos = graphviz_layout(G, prog='dot')
# labels = nx.get_node_attributes(G, 'label')
# #print(pos)
# edge_labels = nx.get_edge_attributes(G, 'label')
# node_sizes = [size for (_, size) in nx.get_node_attributes(G, 'node_size').items()]
# # node_shapes = [shape for (_, shape) in nx.get_node_attributes(G, 'node_shape').items()]
# #node_shapes = nx.get_node_attributes(G, 'node_shape')
#
# G = G.reverse()
# #f = plt.figure()
# # width = 6
# width = 3
# height = 2
# plt.figure(figsize=(width, height))
# pos = {1: (240, 25.0), 2: (260, 25.0), 3: (250, 54.0)}
# # pos = {1: (5, 5), 2: (5, 5.0), 3: (5, 5)}
# nx.draw(G, pos, with_labels=False, arrows=True, node_color='orchid', node_size=node_sizes,font_size=4)
# # nx.draw_networkx_labels(G, pos, with_labels=False, arrows=True, node_color='orchid',
# #                         node_size=node_sizes,font_size=8)
#
# text = nx.draw_networkx_labels(G, pos=pos, font_size=8, labels=labels, rotation='horizontal')
# # for i, t in text.items():
# #     t.set_rotation(90)
#
# # nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, bbox=dict(alpha=0), font_size=7)
# print(pos)
# label_pos = {node: get_label_pos(G, node, pos) for (node, position) in pos.items()}
# labels = {node: get_label(G, node) for (node, position) in pos.items()}
# # print(label_pos)
# # label_pos = {1: get_label_pos(G, 1, pos)}
# text = nx.draw_networkx_labels(G, pos=label_pos, font_size=12, labels=labels, rotation='horizontal')
# for i, t in text.items():
#     if t.get_text()[0:2] in ['Ns', 'sc', 'dc', "$\\", r"$\""]:
#         for k in G.predecessors(i):
#             diff = (label_pos[i][0]-label_pos[k][0], label_pos[i][1]-label_pos[k][1])
#             angle = degrees(atan2(diff[1], diff[0]))
#             print(angle)
#         t.set_rotation('vertical')
#     else:
#         t.set_rotation('vertical')
#         # pass
# #nxpd.draw(G)
# # plt.show()
#
# fig = plt.gcf()
# ax = plt.gca()
# print(ax.get_xlim())
# print(ax.get_ylim())
# ax.set_xlim((235, 265))
# ax.set_ylim((18, 64))
# size = fig.get_size_inches()
# # plt.gcf().su#bplots_adjust(bottom=0.15)
# plt.savefig('Notes/scheme_graphs/distill_graph.pdf')#, bbox_inches='tight')

#
#

#
#
G = nx.DiGraph()



G.add_node(1, node_size=1800, label='ELG')
G.add_node(2, node_size=1800, label='ELG')
G.add_node(3, node_size=1800, label='A-$\\textrm{QR}$')
G.add_node(4, node_size=1800, label='A-$\\textrm{QR}$')
G.add_node(5, node_size=1800, label='A-$\\textrm{QR}$')
G.add_node(6, node_size=1800, label='ELG')
G.add_node(7, node_size=1800, label='$\\textrm{QR}$-B')
G.add_node(8, node_size=1800, label='A-B')

G.add_edge(8,7, label='Swap,\nr = r$_5$')
G.add_edge(7,6, label='ELG,\nr = r$_4$')
G.add_edge(8,5, label='Swap,\nr = r$_5$')
G.add_edge(5,3, label='Distill,\nr = r$_3$')
G.add_edge(5,4, label='Distill,\nr = r$_3$')
G.add_edge(3,1, label='ELG,\nr = r$_1$')
G.add_edge(4,2, label='ELG,\nr = r$_2$')

#save_and_open_pdf_of_scheme_graph(G, 'test11june.png')
pos = graphviz_layout(G, prog='dot')
labels = nx.get_node_attributes(G, 'label')
#print(pos)
edge_labels = nx.get_edge_attributes(G, 'label')
node_sizes = [size for (_, size) in nx.get_node_attributes(G, 'node_size').items()]
# node_shapes = [shape for (_, shape) in nx.get_node_attributes(G, 'node_shape').items()]
#node_shapes = nx.get_node_attributes(G, 'node_shape')

G = G.reverse()
#f = plt.figure()
# width = 6
width = 3
height = 4.3
plt.figure(figsize=(width, height))
# pos = {1: (240, 25.0), 2: (260, 25.0), 3: (250, 54.0)}
# pos = {1: (5, 5), 2: (5, 5.0), 3: (5, 5)}
nx.draw(G, pos, with_labels=False, arrows=True, node_color='orchid', node_size=node_sizes,font_size=4)
# nx.draw_networkx_labels(G, pos, with_labels=False, arrows=True, node_color='orchid',
#                         node_size=node_sizes,font_size=8)

text = nx.draw_networkx_labels(G, pos=pos, font_size=8, labels=labels, rotation='horizontal')
# for i, t in text.items():
#     t.set_rotation(90)

# nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, bbox=dict(alpha=0), font_size=7)
print(pos)
label_pos = {node: get_label_pos(G, node, pos) for (node, position) in pos.items()}
labels = {node: get_label(G, node) for (node, position) in pos.items()}
# print(label_pos)
# label_pos = {1: get_label_pos(G, 1, pos)}
text = nx.draw_networkx_labels(G, pos=label_pos, font_size=12, labels=labels, rotation='horizontal')
for i, t in text.items():
    if t.get_text()[0:2] in ['Ns', 'sc', 'dc', "$\\", r"$\""]:
        for k in G.predecessors(i):
            diff = (label_pos[i][0]-label_pos[k][0], label_pos[i][1]-label_pos[k][1])
            angle = degrees(atan2(diff[1], diff[0]))
            print(angle)
        t.set_rotation('vertical')
    else:
        t.set_rotation('vertical')
        # pass
#nxpd.draw(G)
# plt.show()

fig = plt.gcf()
ax = plt.gca()
print(ax.get_xlim())
print(ax.get_ylim())
ax.set_xlim((34.5, 399.5))
# ax.set_ylim((18, 64))
size = fig.get_size_inches()
# plt.gcf().su#bplots_adjust(bottom=0.15)
# plt.show()
plt.savefig('Notes/scheme_graphs/scheme_graph.pdf')#, bbox_inches='tight')
