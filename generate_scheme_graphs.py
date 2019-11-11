from ete3 import Tree, TreeNode
from plot_functions import generate_pdf_of_scheme_graph, select_data
from plot_functions import save_and_open_pdf_of_scheme_graph
import networkx as nx
from random import randint
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import matplotlib
from matplotlib.figure import figaspect
matplotlib.use('TkAgg')
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


all_schemes = select_data()
max_length_schemes = all_schemes[max(all_schemes)]

# for i, j in max_length_schemes.items():
#     print(i)
#     print(j)

highest_fidelity_scheme = max_length_schemes[max(max_length_schemes, key=lambda x: x[0])]
lowest_fidelity_scheme = max_length_schemes[min(max_length_schemes, key=lambda x: x[0])]

min_key = min(max_length_schemes, key=lambda x: x[0])
max_key = max(max_length_schemes, key=lambda x: x[0])

print(min_key)
print(max_key)

x = []
y = []

highest_fidelity_scheme = max_length_schemes[(46, 0)]
print(highest_fidelity_scheme.fidelity)
print(highest_fidelity_scheme.time)

# fig = plt.figure()
# ax = plt.gca()
# ax.scatter(x, y)
# ax.set_yscale('log')
# plt.show()
G = generate_pdf_of_scheme_graph(highest_fidelity_scheme, filename="test_7th_of_june",
                                 general_distillation=True, node_ID=randint(0, 99999999))



#save_and_open_pdf_of_scheme_graph(G, 'test11june.png')
pos = graphviz_layout(G, prog='dot')
pos = {node: (position[0], position[1]*1) for (node, position) in pos.items()}
print(pos)
labels = nx.get_node_attributes(G, 'label')

edge_labels = nx.get_edge_attributes(G, 'label')
node_sizes = [size for (_, size) in nx.get_node_attributes(G, 'node_size').items()]
node_shapes = [shape for (_, shape) in nx.get_node_attributes(G, 'node_shape').items()]
#print(node_shapes)



G = G.reverse()
# width = 6
width = 11
height = 8
plt.figure(figsize=(width, height))
# nx.draw(G, pos, labels=labels, node_list=G.nodes(), with_labels=labels, arrows=True, node_color='mediumorchid',
#         node_size=node_sizes, font_size=8)

nx.draw(G, pos, with_labels=False, arrows=True, node_color='mediumorchid', node_size=node_sizes,font_size=8)
# nx.draw_networkx_labels(G, pos, with_labels=False, arrows=True, node_color='mediumorchid',
#                         node_size=node_sizes,font_size=8)
text = nx.draw_networkx_labels(G, pos=pos, font_size=8, labels=labels, rotation='vertical')
for i, t in text.items():
    t.set_rotation(90)
pos = {node: (position[0], position[1]) for (node, position) in pos.items()}
# nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, bbox=dict(alpha=0), font_size=7)

label_pos = {node: get_label_pos(G, node, pos) for (node, position) in pos.items()}
labels = {node: get_label(G, node) for (node, position) in pos.items()}
print(label_pos)
text = nx.draw_networkx_labels(G, pos=label_pos, font_size=8, labels=labels, rotation='vertical')
for i, t in text.items():
    if t.get_text()[0:2] in ['Ns', 'sc', 'dc', "$\\", r"$\\"]:
        for k in G.predecessors(i):
            diff = (label_pos[i][0]-label_pos[k][0], label_pos[i][1]-label_pos[k][1])
            angle = degrees(atan2(diff[1], diff[0]))
        t.set_rotation(angle)
    else:
        t.set_rotation('vertical')
        # pass
#nxpd.draw(G)
# plt.show()
fig = plt.gcf()
size = fig.get_size_inches()

plt.savefig('Notes/scheme_graphs/F=096049915_T=00177076_800km_10rptrs_IP_and_MP_sym.pdf', bbox_inches='tight')


G = generate_pdf_of_scheme_graph(lowest_fidelity_scheme, filename="test_7th_of_june",
                                 general_distillation=True, node_ID=randint(0, 99999999))


plt.clf()
#save_and_open_pdf_of_scheme_graph(G, 'test11june.png')
pos = graphviz_layout(G, prog='dot')
labels = nx.get_node_attributes(G, 'label')
print(pos)
edge_labels = nx.get_edge_attributes(G, 'label')
node_sizes = [size for (_, size) in nx.get_node_attributes(G, 'node_size').items()]
node_shapes = [shape for (_, shape) in nx.get_node_attributes(G, 'node_shape').items()]
#node_shapes = nx.get_node_attributes(G, 'node_shape')

G = G.reverse()
#f = plt.figure()
# width = 6
width = 8
height = 8
plt.figure(figsize=(width, height))
nx.draw(G, pos, with_labels=False, arrows=True, node_color='mediumorchid', node_size=node_sizes,font_size=8)
# nx.draw_networkx_labels(G, pos, with_labels=False, arrows=True, node_color='mediumorchid',
#                         node_size=node_sizes,font_size=8)
text = nx.draw_networkx_labels(G, pos=pos, font_size=8, labels=labels, rotation='vertical')
for i, t in text.items():
    t.set_rotation(90)
pos = {node: (position[0], position[1]*1) for (node, position) in pos.items()}
# nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, bbox=dict(alpha=0), font_size=7)

label_pos = {node: get_label_pos(G, node, pos) for (node, position) in pos.items()}
labels = {node: get_label(G, node) for (node, position) in pos.items()}
print(label_pos)
text = nx.draw_networkx_labels(G, pos=label_pos, font_size=8, labels=labels, rotation='vertical')
for i, t in text.items():
    if t.get_text()[0:2] in ['Ns', 'sc', 'dc', "$\\", r"$\""]:
        for k in G.predecessors(i):
            diff = (label_pos[i][0]-label_pos[k][0], label_pos[i][1]-label_pos[k][1])
            angle = degrees(atan2(diff[1], diff[0]))
            print(angle)
        t.set_rotation(angle)
    else:
        t.set_rotation('vertical')
        # pass
#nxpd.draw(G)
# plt.show()
fig = plt.gcf()
size = fig.get_size_inches()

plt.savefig('Notes/scheme_graphs/lowest_800km_10rptrs_IP_and_MP_sym.pdf', bbox_inches='tight')
#plt.clf()


# scheme = lowest_fidelity_scheme

# G = generate_pdf_of_scheme_graph(scheme, filename="test_7th_of_june",
#                                  general_distillation=True, node_ID = randint(0, 99999999))



# #save_and_open_pdf_of_scheme_graph(G, 'test11june.png')
# pos = graphviz_layout(G, prog='dot')
# labels = nx.get_node_attributes(G, 'label')
# print(pos)
# edge_labels = nx.get_edge_attributes(G, 'label')
# node_sizes = [size for (_, size) in nx.get_node_attributes(G, 'node_size').items()]
# node_shapes = [shape for (_, shape) in nx.get_node_attributes(G, 'node_shape').items()]
# #node_shapes = nx.get_node_attributes(G, 'node_shape')
#
# G = G.reverse()
# nx.draw(G, pos, labels=labels, node_list=G.nodes(), with_labels=labels, arrows=True, node_color='mediumorchid',
#         node_size=node_sizes)
# nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)
# #nxpd.draw(G)
# # plt.show()
# plt.savefig('lowest.png')
# plt.clf()