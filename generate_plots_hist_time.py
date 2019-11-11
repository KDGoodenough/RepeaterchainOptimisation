
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

#colors = ['indianred', 'mediumseagreen', 'steelblue', 'gold', 'orchid']
colors = ['mediumseagreen', 'orchid', 'goldenrod', 'cornflowerblue']#, 'orangered', 'mediumblue']


labels = [r'$\varepsilon_\textrm{swap} = \varepsilon_\textrm{distill} = 0.0125$',
            r'$\varepsilon_\textrm{swap} = \varepsilon_\textrm{distill} = 0.025$',
            r'$\varepsilon_\textrm{swap} = \varepsilon_\textrm{distill} = 0.05$',
            r'$\varepsilon_\textrm{swap} = \varepsilon_\textrm{distill} = 0.1$',
            r'$\varepsilon_\textrm{swap} = \varepsilon_\textrm{distill} = 0.2$']


colors = ['mediumseagreen', 'goldenrod', 'orchid','orangered', 'cornflowerblue']
#colors = ['mediumseagreen', 'goldenrod', 'orangered', 'mediumblue']#, 'orchid']#, 'orangered', 'cornflowerblue', 'red']



#w, h = figaspect(1.35)

w, h = figaspect(1.95)

fig, ax = plt.subplots()
#fig, axes = plt.subplots(nrows=4, ncols=1, sharex = True)
#fig.subplots_adjust(hspace.5)

times1 = [98.74074292182922, 108.3717131614685, 132.46123433113098, 173.89092421531677, 240.1660840511322, 296.6190619468689]
times2 = [104.61022305488586, 115.27117395401001, 129.71443605422974, 165.68303799629211, 206.9744701385498, 316.2775752544403]


times1 = [98.74074292182922, 108.3717131614685, 132.46123433113098, 173.89092421531677, 240.1660840511322]
times2 = [104.61022305488586, 115.27117395401001, 129.71443605422974, 165.68303799629211, 206.9744701385498]

#[44.754966020584106, 57.15931272506714, 141.60315895080566, 886.4988131523132, 61792.28430509567, 0

ind = np.arange(5)    # the x locations for the groups
width = 0.35

barlist=ax.bar(ind-width/2-0.025, times1, width = width)
for i, color in enumerate(colors):
   barlist[i].set_color(color)

barlist=ax.bar(ind+width/2+0.025, times2, width = width)
for i, color in enumerate(colors):
    barlist[i].set_color(color)
    barlist[i].set_edgecolor('white')
    barlist[i].set_linewidth(1)
    barlist[i].set_hatch('/')


ax.set_xticklabels(['', 0.0125, 0.025, 0.05, 0.1, 0.2], fontsize = 16)
ax.set_ylabel('Algorithm time (s)', fontsize = 16)
ax.set_xlabel(r'$\varepsilon_{\textrm{distill}},~\varepsilon_{\textrm{swap}}$', fontsize = 16)

plt.yticks(fontsize = 16)
#plt.show()
name = 'Notes/figs/eps_dist_swap_times.pdf'
#name = 'Notes/figs/300 km, 4 rptrs IP eps dist and swap comp1.pdf'
plt.savefig(name, format='pdf', dpi=1000)#, bbox_extra_artists=(leg1,), bbox_inches='tight')