import matplotlib
from matplotlib.figure import figaspect
matplotlib.use('TkAgg')
matplotlib.rc('text', usetex=True)
from plot_functions import *
import warnings
import matplotlib.patches as mpatches

warnings.filterwarnings("error")
warnings.filterwarnings("ignore", message="elementwise comparison failed; returning scalar instead,"
                                          "but in the future will perform elementwise comparison")

labels = ['No node', 'Single node', 'Two nodes']

colors = ['mediumseagreen', 'orchid', 'goldenrod']
linestyles = ['-', '-', '-', '-.', '-.', '-.', ':', ':', ':']*99

patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(labels))]

w, h = figaspect(1.95)

fig, axes = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(w,h))

select_and_generate_plots(colors=colors, labels=labels, linestyles=linestyles, ax=axes[0])
select_and_generate_plots(colors=colors, labels=labels, linestyles=linestyles, ax=axes[1])
select_and_generate_plots(colors=colors, labels=labels, linestyles=linestyles, ax=axes[2])
select_and_generate_plots(colors=colors, labels=labels, linestyles=linestyles, ax=axes[3])

x = plt.xlim()
y = plt.ylim()

plt.subplots_adjust(hspace=.01, bottom=0.01, top=1)

leg1 = axes[0].legend(handles=patches, prop={'size': 9.1}, bbox_to_anchor=(0.5, 1.12),
                      ncol=3, loc='center', handletextpad=0.4)

axes[0].add_artist(leg1)

parameterlines1, = axes[0].plot([-1], [0], '-',  alpha=1, linewidth=2, color='k')
parameterlines2, = axes[0].plot([-1], [0], '-.', alpha=1, linewidth=2, color='k')
parameterlines3, = axes[0].plot([-1], [0], ':',  alpha=1, linewidth=2, color='k')
parameterlines = [parameterlines1, parameterlines2, parameterlines3]

leg2 = axes[0].legend(parameterlines, ['Set 2', 'Set 3', 'Set 4'], bbox_to_anchor=(0.5, 1.42),
                      loc='upper center', prop={'size': 9.2}, ncol=3)

plt.xlim = x
plt.ylim = y

for a in axes:
    a.set_ylim([2e-4, 2e-1])
    a.set_xlim([0.5, 1])

axes[3].set_xlabel('Fidelity', fontsize=12)
plt.ylabel('Generation time (s)', fontsize=12)
plt.setp([a.get_xticklabels() for a in axes[:-1]], visible=False)

for a in axes:
    plt.setp(a.get_xticklabels(), fontsize=13)
    plt.setp(a.get_yticklabels(), fontsize=13)

ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]
ax4 = axes[3]

ax1.get_shared_x_axes().join(ax1, ax2)
ax3.get_shared_x_axes().join(ax1, ax3)

x1, x2, y = 0.025, 0.83, 0.95

ax1.text(x2, y, '(a)', fontsize=12, verticalalignment='top', transform=axes[0].transAxes)
ax2.text(x2, y, '(b)', fontsize=12, verticalalignment='top', transform=axes[1].transAxes)
ax3.text(x2, y, '(c)', fontsize=12, verticalalignment='top', transform=axes[2].transAxes)
ax4.text(x2, y, '(d)', fontsize=12, verticalalignment='top', transform=axes[3].transAxes)

ax1.text(x1, y, '$L$ = 50 km', fontsize=12, verticalalignment='top', transform=axes[0].transAxes)
ax2.text(x1, y, '$L$ = 100 km', fontsize=12, verticalalignment='top', transform=axes[1].transAxes)
ax3.text(x1, y, '$L$ = 150 km', fontsize=12, verticalalignment='top', transform=axes[2].transAxes)
ax4.text(x1, y, '$L$ = 200 km', fontsize=12, verticalalignment='top', transform=axes[3].transAxes)

plt.show()

# plt.savefig(name, format='pdf', dpi=1000)#, bbox_extra_artists=(leg1,), bbox_inches='tight')
