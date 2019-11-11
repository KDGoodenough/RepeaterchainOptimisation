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
    "ignore", message="elementwise comparison failed; returning scalar instead,"
                      " but in the future will perform elementwise comparison")


w, h = figaspect(1.95)

#fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (w,h))
fig, ax = plt.subplots(nrows=1, ncols=1)

generate_param_colour_plot_throughput(ax)
#generate_param_colour_plot_throughput(fig, axes[1])
#generate_param_colour_plot_throughput(fig, axes[2])
#generate_param_colour_plot(0.75, fig, axes[1])
#generate_param_colour_plot(0.8, fig, axes[2])


#axes[2].tick_params(axis='x', bottom=True, top=False)

ax.xaxis.set_tick_params(labeltop=False, labelbottom=True)
ax.set_xlabel('Gate fidelities', fontsize = 18)
#ax.set_ylabel('Success probabilities', fontsize = 18)
#ax.set_xlabel('Number of modes', fontsize = 22)
#ax.set_xlabel('Success probabilities', fontsize = 22)
ax.set_ylabel('Success probabilities', fontsize = 18)
ax.xaxis.set_label_position('bottom')

#ax.set_ylabel('Coherence times', fontsize = 18)
#plt.ylabel('Coherence times', fontsize = 18)


# ax.text(0.88, 0.9, '(a)', fontsize=12,
#         verticalalignment='top', transform = axes[0].transAxes)
#
# ax.text(0.88, 0.9, '(b)', fontsize=12,
#         verticalalignment='top', transform = axes[1].transAxes, color = 'black')
#
# ax.text(0.88, 0.9, '(c)', fontsize=12,
#         verticalalignment='top', transform = axes[2].transAxes, color = 'black')

#name = 'Notes/figs/IP_color_contour, 50 km, 1 rptrs, G=0.98-1 in 20, T=1-100 in 20, eps = 0.01, IP, 1 node, 50km, 8 feb2.pdf'
#name = 'Notes/figs/NP_color_contour, 50 km, 1 rptrs 16 number_of_modes=[1e3,1e4,1e5, 1e6], T=[1e-2, 1e-1, 1e0, 1e1], eps = 0.01, NP.pdf'
#name = 'Notes/figs/IP_contour_SKR_sixstate_changing_P.pdf'
#name = 'Notes/figs/MP_contour_15_km_SKR_sixstate_modes+Tmarch15.pdf'
ax.set_yticklabels(['', '1e-2', '5.6e-2', '4.4e-1', '2.5e0', '1e1'])
ax.set_xticklabels(['', '1e3', '5.6e3', '4.4e4', '2.5e5', '1e6'])
ax.set_yticklabels(['', '1e-2', '5.6e-2', '4.4e-1', '2.5e0', '1e1'], fontsize = 15)
ax.set_xticklabels(['', '1e4', '5.6e4', '4.4e5', '2.5e6', '1e7'], fontsize = 15)
#ax.set_yticklabels(['', r'$1.0\times 10^{-2}$', r'$3.7\times 10^{-2}$', r'$5.7\times 10^{-2}$', r'$8.1\times 10^{-2}$'], fontsize = 15)
#ax.set_xticklabels(['', r'$0.900$', r'$0.913$', r'$0.926$', r'$0.939$'], fontsize = 15)

# ax.set_yticklabels(['', '0.8', '0.85', '0.9', '0.95'])
# ax.set_xticklabels(['', '0.980', '0.985', '0.990', '0.995'])

ax.set_xticklabels(['', 0.98, 0.9822, 0.9844, 0.9866, 0.9888])
ax.set_yticklabels(['', 1, 1.88, 2.77, 3.66, 4.55])
ax.set_xticklabels(['', 0.98, 0.985, 0.99, 0.995, 1])
ax.set_yticklabels(['', 0.8, 0.85, 0.9, 0.95, 1])
ax.set_yticklabels(['', 1, 3.37, 5.74, 8.11, 10], fontsize = 16)
ax.set_yticklabels(['', 0.8, 0.853, 0.905, 0.958, 1], fontsize = 16)
# ax.set_xticklabels(['', '1e3', '1e4', '1e5', '1e6'])
ax.set_xticklabels(['', 0.98, 0.985, 0.991, 0.996], fontsize = 16)


plt.rc('text', usetex=True)
# ax.set_yticklabels(['', 1, 3.37, 5.74, 8.11])
# ax.set_yticklabels(['', 0.8, 0.85, 0.91, 0.96])
# ax.set_xticklabels(['', 0.980, 0.985, 0.991, 0.996, 1])
print(fig)
plt.gcf().subplots_adjust(left=0.4)
# plt.show()
#plt.tight_layout()
#plt.show()
#name = 'Notes/figs/MP_contour_SKR_sixstate_coh+modes+Tmarch15.pdf'
# name = "Notes/figs/IP_contour_SKR_sixstate_21_march.pdf"
name = "Notes/figs/IP_contour_SKR_sixstate_changing_P_21_march.pdf"
plt.savefig(name, format='pdf', dpi=1000)#, bbox_extra_artists=(leg1,), bbox_inches='tight')