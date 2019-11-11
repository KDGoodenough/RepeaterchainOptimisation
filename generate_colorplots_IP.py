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


w, h = figaspect(1.95)

fig, axes = plt.subplots(nrows=3, ncols=1, sharex = True, figsize = (w,h))

generate_param_colour_plot(0.7, fig, axes[0], scale =30)
generate_param_colour_plot(0.8, fig, axes[1], scale = 32)
generate_param_colour_plot(0.9, fig, axes[2], scale = 38)
#generate_param_colour_plot(0.75, fig, axes[1])
#generate_param_colour_plot(0.8, fig, axes[2])


#axes[2].tick_params(axis='x', bottom=True, top=False)

axes[2].xaxis.set_tick_params(labeltop=False, labelbottom=True)
axes[1].set_ylabel('Coherence times', fontsize = 12)
axes[2].set_xlabel('Success probabilities', fontsize = 12)
#axes[2].set_xlabel('Number of modes', fontsize = 18)
axes[2].xaxis.set_label_position('bottom')

#axes[1].set_xlabel('Coherence time', fontsize = 18)
#plt.ylabel('Coherence times', fontsize = 18)

axes[0].text(0.03, 0.975, '(a)', fontsize=12,
        verticalalignment='top', transform = axes[0].transAxes, color = 'black')

axes[1].text(0.03, 0.975, '(b)', fontsize=12,
        verticalalignment='top', transform = axes[1].transAxes, color = 'black')

axes[2].text(0.03, 0.975, '(c)', fontsize=12,
        verticalalignment='top', transform = axes[2].transAxes, color = 'black')

#name = 'Notes/figs/IP_color_contour, 50 km, 1 rptrs, G=0.98-1 in 20, T=1-100 in 20, eps = 0.01, IP, 1 node, 50km, 8 feb2.pdf'
#name = 'Notes/figs/NP_color_contour, 50 km, 1 rptrs 16 number_of_modes=[1e3,1e4,1e5, 1e6], T=[1e-2, 1e-1, 1e0, 1e1], eps = 0.01, NP.pdf'
name = 'Notes/figs/15km-MP_contour_fidelity2_march_21_normed.pdf'
#name = "Notes/figs/IP_color_contour, 50 km, 1 rptrs, G=0.98-1 in 20, P=0.8-1 in 20, eps = 0.01, IP, 1 node, 50km, 19 march.pdf"
#name = "Notes/figs/IP_color_contour_19_march.pdf"
axes[2].set_xticklabels(['', 0.9, 0.913, 0.926, 0.939])
print(fig)
print(axes)
plt.gcf().subplots_adjust(left=0.4)
#plt.show()
#plt.tight_layout()
# plt.show()
plt.savefig(name, format='pdf', dpi=1000)#, bbox_extra_artists=(leg1,), bbox_inches='tight')