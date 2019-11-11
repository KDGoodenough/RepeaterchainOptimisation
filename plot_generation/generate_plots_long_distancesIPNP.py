import matplotlib
from matplotlib.figure import figaspect
matplotlib.use('TkAgg')
from matplotlib import pyplot
matplotlib.rc('text', usetex = True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from plot_functions import *
import time
from itertools import product
import warnings
import plot_functions

import matplotlib.patches as mpatches
warnings.filterwarnings("error")
warnings.filterwarnings(
    "ignore", message="elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison")







#for i in range(2, 10):
#    print(i, list(create_k_range(True, i)))
#    print(i, calc_g_h(i))
#    print(i, list(create_k_range(True, i)))
#    print('--')
#print(time2-time1)
#print(time3-time2)
# fidelities, times = gather_schemes_from_files()

#import matplotlib.pyplot as plt

#time = [39.67356085777283, 51.28681683540344, 371.49295806884766, 972.099543094635, 6399.4047498703]
#epsilons = [0.1, 0.05, 0.01, 0.005, 0.001]

colors = ['indianred', 'mediumseagreen', 'steelblue', 'gold', 'orchid']


# plt.step(fidelities[0], times[0], '-', where='pre',
#              linewidth=4, ms=18, color=colors[0], alpha=0.8, label = r'No nodes')
# plt.step(fidelities[1], times[1], '-', where='pre',
#               linewidth=4, ms=18, color=colors[1], alpha=0.8, label = r'Single node')
# plt.step(fidelities[2], times[2], '-', where='pre',
#                linewidth=4, ms=18, color=colors[2], alpha=0.8, label = r'Two nodes')
# plt.step(fidelities[3], times[3], '-', where='pre',
#               linewidth=4, ms=18, color=colors[3], alpha=0.8, label = r'$\varepsilon_F = 0.005$')
# plt.step(fidelities[4], times[4], '-', where='pre',
#               linewidth=4, ms=18, color=colors[4], alpha=0.8, label = r'$\varepsilon_F = 0.001$')
# # plt.plot(epsilons, time, linewidth = 2, linestyle = ':', marker="h", markersize = 8)
# # plt.plot(epsilons, time, markersize = 5, marker= "h")
# plt.xlim(0.5, 0.853)
# plt.ylim(0, 0.123)
# plt.xlabel(r'$F$', fontsize=32)
# plt.ylabel('Time (s)', fontsize=32)
# plt.title(r'Distance of 50 km', fontsize=32)
# plt.xticks(fontsize=30)
# plt.yticks(fontsize=30)
# # #plt.legend((r'$\varepsilon_F = 0.1',r'$\varepsilon_F = 0.1', r'$\varepsilon_F = 0.1', r'$\varepsilon_F = 0.1', r'$\varepsilon_F = 0.1'))
# plt.legend(loc='best')
# plt.legend(prop={'size':25})
# plt.show()




#plt.step(fidelities[0], times[0], '-', where='pre',
#              linewidth=4, ms=18, color=colors[0], alpha=0.8, label = r'$\varepsilon_d = 1/0.1, \varepsilon_s = 1/0.1,$')
# plt.step(fidelities[1], times[1], '-', where='pre',
#               linewidth=4, ms=18, color=colors[1], alpha=0.8, label = r'$\varepsilon_F = 0.05$')
# plt.step(fidelities[2], times[2], '-', where='pre',
#               linewidth=4, ms=18, color=colors[2], alpha=0.8, label = r'$\varepsilon_F = 0.01$')
# plt.step(fidelities[3], times[3], '-', where='pre',
#               linewidth=4, ms=18, color=colors[3], alpha=0.8, label = r'$\varepsilon_F = 0.005$')
# plt.step(fidelities[4], times[4], '-', where='pre',
#               linewidth=4, ms=18, color=colors[4], alpha=0.8, label = r'$\varepsilon_F = 0.001$')
# # plt.plot(epsilons, time, linewidth = 2, linestyle = ':', marker="h", markersize = 8)
# # plt.plot(epsilons, time, markersize = 5, marker= "h")
# plt.xlabel(r'$F$', fontsize=32)
# plt.ylabel('Time (s)', fontsize=32)
# plt.title(r'Optimised schemes as a function of ' r'$\varepsilon_F$', fontsize=32)
# plt.xticks(fontsize=30)
# plt.yticks(fontsize=30)
# #plt.legend((r'$\varepsilon_F = 0.1',r'$\varepsilon_F = 0.1', r'$\varepsilon_F = 0.1', r'$\varepsilon_F = 0.1', r'$\varepsilon_F = 0.1'))
# plt.legend(loc='best')
# plt.legend(prop={'size':20})
# plt.show()







#plt.step(fidelities[0], times[0], '-', where='pre',
#              linewidth=4, ms=18, color=colors[0], alpha=0.8, label = r'$\varepsilon_d = 1/0.1, \varepsilon_s = 1/0.1,$')
# plt.step(fidelities[1], times[1], '-', where='pre',
#               linewidth=4, ms=18, color=colors[1], alpha=0.8, label = r'$\varepsilon_F = 0.05$')
# plt.step(fidelities[2], times[2], '-', where='pre',
#               linewidth=4, ms=18, color=colors[2], alpha=0.8, label = r'$\varepsilon_F = 0.01$')
# plt.step(fidelities[3], times[3], '-', where='pre',
#               linewidth=4, ms=18, color=colors[3], alpha=0.8, label = r'$\varepsilon_F = 0.005$')
# plt.step(fidelities[4], times[4], '-', where='pre',
#               linewidth=4, ms=18, color=colors[4], alpha=0.8, label = r'$\varepsilon_F = 0.001$')
# # plt.plot(epsilons, time, linewidth = 2, linestyle = ':', marker="h", markersize = 8)
# # plt.plot(epsilons, time, markersize = 5, marker= "h")
# plt.xlabel(r'$F$', fontsize=32)
# plt.ylabel('Time (s)', fontsize=32)
# plt.title(r'Optimised schemes as a function of ' r'$\varepsilon_F$', fontsize=32)
# plt.xticks(fontsize=30)
# plt.yticks(fontsize=30)
# #plt.legend((r'$\varepsilon_F = 0.1',r'$\varepsilon_F = 0.1', r'$\varepsilon_F = 0.1', r'$\varepsilon_F = 0.1', r'$\varepsilon_F = 0.1'))
# plt.legend(loc='best')
# plt.legend(prop={'size':20})
# plt.show()



labels = ['Single node', 'Three nodes', 'Five nodes']
labels = ['No node', 'Single node', 'Two nodes']
labels = ['No node', 'Single node, $N = 0$', 'Single node, $N = 1$', 'Single node, $N = 2$', 'Single node, $N = 3$']
labels = ['No node, 20 km', 'Single node, 40 km', 'Two nodes, 60 km', 'Three nodes, 80 km', 'Four nodes, 100 km', 'Five nodes, 120 km']
labels = ['No node', 'Single node', 'Two nodes']
labels = ['Fifteen nodes', 'Twenty nodes', 'Thirty nodes', 'Forty nodes']
#labels = ['No node', 'Single node']
#labels = ['No node', 'Single node', 'Two nodes', 'Three nodes', 'Four nodes', 'Five nodes']
#labels = ['No node', 'Two nodes', 'Four nodes', 'Six nodes']
#labels = ['', r'$\varepsilon_F = 1\mathrm{e}{-1}$,' + '\n' + r'$ \varepsilon_p = 2\mathrm{e}{-1}$', r'$\varepsilon_F = 5\mathrm{e}{-2}$,' + '\n' + r'$\varepsilon_p = 1\mathrm{e}{-1}$',
#          r'$\varepsilon_F = 1\mathrm{e}{-2}$,' + '\n' + r'$ \varepsilon_p = 2\mathrm{e}{-2}$', r'$\varepsilon_F = 5\mathrm{e}{-3}$,' + '\n' + r'$\varepsilon_p = 1\mathrm{e}{-2}$', r'$\varepsilon_F = 1\mathrm{e}{-3}$,' + '\n' + r'$\varepsilon_p = 2\mathrm{e}{-3}$']
colors = ['mediumseagreen', 'orchid', 'goldenrod', 'cornflowerblue', 'orangered']#, 'cornflowerblue', 'orangered', 'mediumblue']#, 'orchid']#, 'orangered', 'cornflowerblue', 'red']
colors = ['orchid', 'goldenrod', 'cornflowerblue', 'orangered']#, 'cornflowerblue', 'orangered', 'mediumblue']#, 'orchid']#, 'orangered', 'cornflowerblue', 'red']
#colors = ['mediumseagreen', 'goldenrod', 'orangered', 'mediumblue']#, 'orchid']#, 'orangered', 'cornflowerblue', 'red']
linestyles = ['-', '-', '-.', '-.', ':', ':']
linestyles = ['-', '-', '-', '-.', '-.', '-.', ':', ':', ':']*99999
linestyles = ['-', ':']*99999
#linestyles = [':']*9999
linestyles = None

patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(labels))]
#
#
#select_and_generate_plots(colors = colors, labels = labels)


#w, h = figaspect(1.35)

w, h = figaspect(1.95)

fig, axes = plt.subplots()
#fig, axes = plt.subplots(nrows=4, ncols=1, sharex = True)
#fig.subplots_adjust(hspace.5)


select_and_generate_plots(colors = colors, labels = labels, linestyles = linestyles, ax=axes)

#select_and_generate_plots(colors = colors, labels = labels, linestyles = linestyles, ax=axes[1])

#select_and_generate_plots(colors = colors, labels = labels, linestyles = linestyles, ax=axes[2])

#select_and_generate_plots(colors = colors, labels = labels, linestyles = linestyles, ax=axes[3])

#legend1 = pyplot.legend(['-', '--', ':'], ["algo1", "algo2", "algo3"], loc=1)

#fig = plt.gcf()
#ax = plt.gca()

# Add second legend for the maxes and mins.
# leg1 will be removed from figure
x = plt.xlim()
y = plt.ylim()


plt.subplots_adjust(hspace=.01, bottom=0.01, top = 1)
#plt.subplots_adjust(hspace=.01, bottom=-0.2)

#plt.gcf().canvas.draw()
#p = leg1.get_window_extent()
#print(p)

# parameterlines1, = axes[0].plot([-1], [0], '-',  alpha=1, linewidth = 2, color = 'k')
# parameterlines2, = axes[0].plot([-1], [0], '-.', alpha=1, linewidth = 2, color = 'k')
# parameterlines3, = axes[0].plot([-1], [0], ':',  alpha=1, linewidth = 2, color = 'k')





# parameterlines = [parameterlines1, parameterlines2, parameterlines3]

# leg1 = axes[0].legend(parameterlines, ['Set 3', 'Set 4'], bbox_to_anchor=(0.5, 1.2), loc=3,
#           ncol=3, mode="expand", borderaxespad=0., prop={'size': 9.5}, handlelength=4, handletextpad=0.4)


# leg1 = axes[0].legend(parameterlines, ['Set 2', 'Set 3', 'Set 4'], bbox_to_anchor=(0.5, 1.28), loc='center',
#           ncol=3, borderaxespad=0., prop={'size': 9.2}, handlelength=3.5, handletextpad=1.0)

#leg1 = axes[0].legend(parameterlines, ['Set 2','Set 3', 'Set 4'], bbox_to_anchor=(0.45, 1.2), loc=3,
#           ncol=3, mode="expand", borderaxespad=0., prop={'size': 14.5}, handlelength=4, handletextpad=0.4)




#leg1 = axes[0].legend(parameterlines, ['Set 2','Set 3', 'Set 4'], bbox_to_anchor=(0.45, 1.34), loc='center',
#           ncol=3, prop={'size': 9.5}, handlelength=3.61, handletextpad=0.95)



#axes[0].legend(handles=patches, prop={'size': 9.2}, bbox_to_anchor=(0.5, 1.10),
 #       ncol=3, loc = 'center', handletextpad=0.4)
axes.legend(handles=patches, prop={'size': 10.8}, bbox_to_anchor=(0.5, 1.05),
        ncol=4,loc = 'center', handletextpad=0.35)

#axes[0].legend(handles=patches, prop={'size': 11.5}, bbox_to_anchor=(0.5, 1.10),
 #        ncol=3, loc = 'center', handletextpad=0.4)



name = 'Notes/figs/test.pdf'

# Manually add the first legend back
# axes[0].add_artist(leg1)

#leg2 = ax.legend(parameterlines,['Set 2','Set 3'], bbox_to_anchor=(0.5, -0.05), loc='best', prop={'size': 14})
#leg2 = plt.legend(handles=patches, prop={'size': 16}, loc='best', ncol=1)


# Put a legend to the right of the current axis


#plt.legend(handles=patches, bbox_to_anchor=(1.05, 2.75), loc=2, borderaxespad=0)



plt.xlim = x
plt.ylim = y

# fig, ax = plt.subplots()
#
# barlist=ax.bar([1,2,3,4,5, 6], [44.754966020584106, 57.15931272506714, 141.60315895080566, 886.4988131523132, 61792.28430509567, 0])
# for i, color in enumerate(colors):
#     barlist[i].set_color(color)
#
#
#
# ax.set_xticklabels(tuple(labels))


#
#

# for a in axes:
#     a.set_xlim([0.5, 1])
    #a.set_ylim([6.5e-5, 0.005])
    #   a.set_ylim([1e-4, 10])

#for a in axes:
    #a.xlabel('Fidelity', fontsize=18)
    #a.ylabel('Generation time (s)', fontsize=18)
    #plt.ylabel('Algorithm time (s)', fontsize=18)
    #plt.xlim((0.25,5.75))
    #plt.setp(a.get_yticklabels(), fontsize=14)
    #plt.setp(a.get_xticklabels(), fontsize=14)
    #a.set_xticklabels([])
    #a.set_yticklabels([])
    #plt.ylim((3*10**(-4), 4*10**(-3)))

axes.set_xlabel('Fidelity', fontsize=18)
#
#axes[2].set_ylabel('Generation time (s)', fontsize=12)
plt.ylabel('Generation time (s)', fontsize=18)
#axes[2].yaxis.set_label_coords(-0.08, 1.5)
#axes[0].set_xticklabels(fontsize = 0)
#axes[1].set_xticklabels(fontsize = 0)
#plt.setp(axes[2].get_xticklabels(), fontsize = 14)

#a.xlabel('Fidelity', fontsize=18)
    #a.ylabel('Generation time (s)', fontsize=18)


# plt.xlabel('Fidelity', fontsize=18)
# plt.ylabel('Generation time (s)', fontsize=18)
# #plt.ylabel('Algorithm time (s)', fontsize=18)
# #plt.xlim((0.25,5.75))
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=18)
# #plt.ylim((3*10**(-4), 4*10**(-3)))

#legend1 = pyplot.legend(handles=patches, loc=4)
#legend2 = pyplot.legend(handles=patches, loc=3)

#plt.add_artist(legend1)
#plt.add_artist(legend2)


# plt.setp([a.get_xticklabels() for a in axes[:-1]], visible=False)


#for a in axes:
plt.setp(axes.get_xticklabels(), fontsize=13)
plt.setp(axes.get_yticklabels(), fontsize=13)
    #a.set_ylim((0.00015, 0.3))



#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width, box.height * 0.8])

# Put a legend to the right of the current axis
#ax.legend(handles=patches, prop={'size': 16.0}, bbox_to_anchor=(0.5, 1.35),
#         ncol=3, loc = 'center')






#axes[1].legend(handles=patches, prop={'size': 16}, loc='best', ncol=1)

#fig.subplots_adjust(hspace=0.01, wspace=0.01, bottom=-0.5, top=1.5)
#fig.subplots_adjust(hspace=-1.0, wspace=0.01, bottom=-0.5, top=1.5)
#plt.subplots_adjust(hspace=-1.0, wspace=0.01, bottom=-0.5, top=1.5)

#plt.subplots_adjust(hspace=0.001, wspace=0.001, bottom=-0.75, top=1.75)



# ax1.text(0.95, 1e-1, '(a)', fontsize=12,
#         verticalalignment='top')
#
# ax2.text(0.95, 1e-1, '(b)', fontsize=12,
#         verticalalignment='top')
#
# ax3.text(0.95, 1e-1, '(c)', fontsize=12,
#         verticalalignment='top')
#
# axes[3].text(0.51, 1e-1, '$L$ = 50 km', fontsize=12,
#         verticalalignment='top')
#
#
# ax1.text(0.51, 1e-1, '$L$ = 50 km', fontsize=12,
#         verticalalignment='top')
#
# ax2.text(0.51, 1e-1, '$L$ = 100 km', fontsize=12,
#         verticalalignment='top')
#
# ax3.text(0.51, 1e-1, '$L$ = 150 km', fontsize=12,
#         verticalalignment='top')
#
# axes[3].text(0.51, 1e-1, '$L$ = 200 km', fontsize=12,
#         verticalalignment='top')



# ax1.text(0.95, 1.5e0, '(a)', fontsize=12,
#         verticalalignment='top')
#
# ax2.text(0.95, 1.5e0, '(b)', fontsize=12,
#         verticalalignment='top')
#
# ax3.text(0.95, 1.5e0, '(c)', fontsize=12,
#         verticalalignment='top')
#
# axes[3].text(0.95, 1.5e0, '(d)', fontsize=12,
#         verticalalignment='top')
#
#
# ax1.text(0.51, 1.5e0, '$L$ = 50 km', fontsize=12,
#         verticalalignment='top')
#
# ax2.text(0.51, 1.5e0, '$L$ = 100 km', fontsize=12,
#         verticalalignment='top')
#
# ax3.text(0.51, 1.5e0, '$L$ = 150 km', fontsize=12,
#         verticalalignment='top')
#
# axes[3].text(0.51, 1.5e0, '$L$ = 200 km', fontsize=12,
#         verticalalignment='top')





x1, x2, y = 0.025, 0.83, 0.95

# ax1.text(x1, y, '$L$ = 50 km', fontsize=12,
#         verticalalignment='top', transform = axes[0].transAxes)
#
# ax2.text(x1, y, '$L$ = 100 km', fontsize=12,
#         verticalalignment='top', transform = axes[1].transAxes)
#
# ax3.text(x1, y, '$L$ = 150 km', fontsize=12,
#         verticalalignment='top', transform = axes[2].transAxes)
#
# ax4.text(x1, y, '$L$ = 200 km', fontsize=12,
#         verticalalignment='top', transform = axes[3].transAxes)





#plt.show()
#plt.savefig('Notes/figs/NVpset2-012rptrs-50km.pdf', format='pdf', dpi=1000m, bbox_extra_artists=(lgd,), bbox_inches='tight')
#plt.savefig('Notes/figs/NVpset2-3-4-012rptrs-200km.pdf', format='pdf', dpi=1000, bbox_extra_artists=(lgd,), bbox_inches='tight')
#plt.savefig('Notes/figs/NVpset4-0246810-200km.pdf', format='pdf', dpi=1000, bbox_extra_artists=(lgd,), bbox_inches='tight')
#plt.savefig('Notes/figs/test.pdf', format='pdf', dpi=1000, bbox_extra_artists=(leg2,), bbox_inches='tight')
#plt.show()

#plt.tight_layout()
#plt.subplot_tool()
#plt.show()
#name = 'Notes/figs/NP-short_distances.pdf'

plt.show()
name = 'Notes/figs/NP-4000km_12_march.pdf'
#plt.savefig(name, format='pdf', dpi=1000)#, bbox_extra_artists=(leg1,), bbox_inches='tight')