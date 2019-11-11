import matplotlib
from matplotlib.figure import figaspect
matplotlib.rc('text', usetex=True)
from plot_functions import *
import warnings
import matplotlib.patches as mpatches
warnings.filterwarnings("error")
warnings.filterwarnings("ignore", message="elementwise comparison failed; returning scalar instead, but in the future "
                                          "will perform elementwise comparison")


labels = ['No node', 'Single node', 'Two nodes']
colors = ['mediumseagreen', 'mediumseagreen', 'orchid', 'orchid', 'goldenrod', 'goldenrod']

linestyles = ['-', ':', '-', ':', '-', ':']
patchcolors = ['mediumseagreen', 'orchid', 'goldenrod']

patches = [mpatches.Patch(color=patchcolors[i], label=labels[i]) for i in range(len(labels))]

w, h = figaspect(1.95)
fig, axes = plt.subplots(nrows=1, ncols=1)

select_and_generate_plots(colors=colors, labels=labels, linestyles=linestyles, ax=axes)
plt.subplots_adjust(hspace=.01, bottom=0.01, top=1)

parameterlines1, = axes.plot([-1], [0], ':',  alpha=1, linewidth=2, color='k')
parameterlines2, = axes.plot([-1], [0], '-', alpha=1, linewidth=2, color='k')
parameterlines = [parameterlines1, parameterlines2]

leg1 = axes.legend(parameterlines, ['With double-click', 'Without double-click'], bbox_to_anchor=(0.5, 1.11), loc=3,
                   ncol=2, mode="expand", borderaxespad=0.7, prop={'size': 11.9}, handlelength=3.5, handletextpad=None,
                   columnspacing=1.4, borderpad=0, frameon=False)

axes.legend(handles=patches, prop={'size': 12.9}, bbox_to_anchor=(0.5, 1.07), ncol=3, loc='center', handletextpad=0.8)

axes.add_artist(leg1)
axes.set_xlim([0.5, 1])
axes.set_xlabel('Fidelity', fontsize=16)

plt.ylabel('Generation time (s)', fontsize=16)

plt.setp(axes.get_xticklabels(), fontsize=15)
plt.setp(axes.get_yticklabels(), fontsize=15)

name = 'IP-50_012-rptrs.pdf'
plt.savefig(name, format='pdf', dpi=1000)
