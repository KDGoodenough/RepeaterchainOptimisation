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


colors = ['indianred', 'mediumseagreen', 'steelblue', 'gold', 'orchid']



labels = ['Single node', 'Three nodes', 'Five nodes']
labels = ['No node', 'Single node', 'Two nodes']
labels = ['No node', 'Single node, $N = 0$', 'Single node, $N = 1$', 'Single node, $N = 2$', 'Single node, $N = 3$']
labels = ['No node, 20 km', 'Single node, 40 km', 'Two nodes, 60 km', 'Three nodes, 80 km', 'Four nodes, 100 km', 'Five nodes, 120 km']
labels = ['200 km', '400 km', '600 km']
#labels = ['No node', 'Single node']
#labels = ['No node', 'Single node', 'Two nodes', 'Three nodes', 'Four nodes', 'Five nodes']
#labels = ['No node', 'Two nodes', 'Four nodes', 'Six nodes']
#labels = ['', r'$\varepsilon_F = 1\mathrm{e}{-1}$,' + '\n' + r'$ \varepsilon_p = 2\mathrm{e}{-1}$', r'$\varepsilon_F = 5\mathrm{e}{-2}$,' + '\n' + r'$\varepsilon_p = 1\mathrm{e}{-1}$',
#          r'$\varepsilon_F = 1\mathrm{e}{-2}$,' + '\n' + r'$ \varepsilon_p = 2\mathrm{e}{-2}$', r'$\varepsilon_F = 5\mathrm{e}{-3}$,' + '\n' + r'$\varepsilon_p = 1\mathrm{e}{-2}$', r'$\varepsilon_F = 1\mathrm{e}{-3}$,' + '\n' + r'$\varepsilon_p = 2\mathrm{e}{-3}$']
colors = ['mediumseagreen', 'orchid', 'goldenrod', 'cornflowerblue']#, 'orangered', 'mediumblue']#, 'orchid']#, 'orangered', 'cornflowerblue', 'red']
#colors = ['mediumseagreen', 'goldenrod', 'orangered', 'mediumblue']#, 'orchid']#, 'orangered', 'cornflowerblue', 'red']
linestyles = ['-', '-', '-.', '-.', ':', ':']
linestyles = ['-', '-', '-', '-.', '-.', '-.', ':', ':', ':']
#linestyles = [':']*9999
linestyles = None

#linestyles = ['-']*6 + [':']*6

patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(labels))]
patch2  = [mpatches.Patch(color='w', label='')]

fig, axes = plt.subplots(nrows=2, ncols=3, sharex = True, sharey = True, figsize = (7,5))
plt.subplots_adjust(hspace=-0.5, wspace = -0.5, bottom=0.0, top = 1)

#select_and_generate_plots(colors = colors, labels = labels, linestyles = linestyles, ax=axes[0])
select_and_generate_plots(colors = colors, labels = labels, linestyles = linestyles, ax=axes[0][0])
select_and_generate_plots(colors = colors, labels = labels, linestyles = linestyles, ax=axes[0][1])
select_and_generate_plots(colors = colors, labels = labels, linestyles = linestyles, ax=axes[0][2])
select_and_generate_plots(colors = colors, labels = labels, linestyles = linestyles, ax=axes[1][0])
select_and_generate_plots(colors = colors, labels = labels, linestyles = linestyles, ax=axes[1][1])
select_and_generate_plots(colors = colors, labels = labels, linestyles = linestyles, ax=axes[1][2])



axes[0][1].legend(handles=patch2
                  , prop={'size': 0.5}, bbox_to_anchor=(0.6, 1.085),
         ncol=1, loc = 'center', handletextpad=0.4, framealpha=0)

lgd = fig.legend(handles=patches, prop={'size': 12.5}, bbox_to_anchor=(0.54, 0.955),
         ncol=4, loc = 'center', handletextpad=0.4)


#axes[0][1].get_legend().remove()
#axes[0][1].add_artist(lgd)
name = 'Notes/figs/parameterexplorationlongdistanceMP10rptrs.pdf'

# Manually add the first legend back
# axes[0].add_artist(leg1)



#axes[3].set_xlabel('Fidelity', fontsize=12)
#
#axes[2].set_ylabel('Generation time (s)', fontsize=12)
axes[1][0].set_ylabel('Generation time (s)', fontsize=18)
axes[1][1].set_xlabel('Fidelity', fontsize=18)

#plt.setp([a.get_xticklabels() for a in axes[:-1]], visible=False)


#for a in axes:
    # plt.setp(a.get_xticklabels(), fontsize=13)
    # plt.setp(a.get_yticklabels(), fontsize=13)
    #a.set_ylim((1.2e-3, 3e0))


ax = axes[0]

ax1 = axes[0]


x1, x2, y = 0.025, 0.915, 0.95

# ax1.text(x2+0.01, y, '(a)', fontsize=12,
#         verticalalignment='top', transform = axes[0].transAxes)
#
# ax2.text(x2-0.15, y, '(b)', fontsize=12,
#         verticalalignment='top', transform = axes[1].transAxes)
#
# ax3.text(x2, y, '(c)', fontsize=12,
#         verticalalignment='top', transform = axes[2].transAxes)
#
# ax4.text(x2, y, '(d)', fontsize=12,
#         verticalalignment='top', transform = axes[3].transAxes)
#
#
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
#
#
#
#


#plt.tight_layout()
#plt.subplot_tool()
#plt.show()
#plt.subplots_adjust(top = 1.5)
plt.savefig(name, format='pdf', dpi=1000)