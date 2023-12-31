"""
* Plot the dynaical-state outputs of the model as a function of a range of amplitude b. Create the figure 8 of the paper.
* The configuration is defined within the "set_up" function in the "set_up.py" script.
* First you need to create the outputs from the "outputs_dynamicalstate_sensitivitytest.py" script.
"""

import matplotlib.pyplot as plt
import numpy as np
import scienceplots
from matplotlib.gridspec import GridSpec
from set_up import set_up_dynamicalstate, set_up_steadystate
import sys
sys.path.append('../')

plt.style.use(['science','no-latex'])
plt.close('all')

# Configuration
dt, end_time, time, Psupply, Psupply_cst, Psupply_senstest, min_param, max_param,l_param = set_up_steadystate()
dt, Psupply_ini, b, nbpulse, end_time_flux, time_flux, nb_time_flux = set_up_dynamicalstate()

## Color map configuration
colors = ['darkgreen', '#ffffff', 'chartreuse']
n_bins = 15
cmap_name = 'custom_green'
cm = plt.cm.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

colors2 = ['#ffffff', 'gray', 'black']
n_bins = 15
cmap_name2 = 'custom_grays'
cm2 = plt.cm.colors.LinearSegmentedColormap.from_list(cmap_name2, colors2, N=n_bins)


plt.rc('font', size=7)

fig = plt.figure(figsize=(10, 7))
gs = GridSpec(3, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1, 1]) 

# Load data from the outputs_dynamicalstate.py code
data1 = np.loadtxt('../outputs/dynamicalstate_pulse1_b0.08_d90.txt')
time_flux1, P1_1, P2_1, Z_1, PO4_1, Psupply1 = data1.T
data2 = np.loadtxt('../outputs/dynamicalstate_pulse2_b0.08_d90.txt')
time_flux2, P1_2, P2_2, Z_2, PO4_2, Psupply2 = data2.T
data3 = np.loadtxt('../outputs/dynamicalstate_pulse3_b0.08_d90.txt')
time_flux3, P1_3, P2_3, Z_3, PO4_3, Psupply3 = data3.T

# Loop to create the subplot

# First line
for i, data in enumerate([(time_flux1, P1_1, P2_1, Z_1, PO4_1, Psupply1),
                           (time_flux2, P1_2, P2_2, Z_2, PO4_2, Psupply2),
                           (time_flux3, P1_3, P2_3, Z_3, PO4_3, Psupply3)]):
    ax = plt.subplot(gs[0, i])
    time_flux, P1, P2, Z, PO4, Psupply = data
    ax.plot(time_flux, Psupply, label=r'$P_{supply}$', color="darkgray")
    ax.plot(time_flux, PO4, label=r'$PO_4$', color="magenta")
    ax.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
    ax.plot(time_flux, P2, label=r'$P_2$', color="green")
    ax.plot(time_flux, Z, label=r'$Z$', color="aqua")
    ax.set_xlabel('Time [d]', fontsize=10)
    ax.set_ylabel('Masses [mmolC m$^{-3}$]', fontsize=10)
    ax.set_ylim(0, 0.65)
    ax.set_xlim(0, 90)
    if i == 0:
        ax.legend(frameon=True, loc='upper left')  # Légende uniquement pour la première figure
        ax.set_title('"One pulse"', fontweight='bold', fontsize=12)
    if i == 1:
        ax.set_title('"Two pulses"', fontweight='bold', fontsize=12)
    if i == 2:
        ax.set_title('"Three pulses"', fontweight='bold', fontsize=12)

# Second line (data from the outputs_dynamicalstate_sensitivitytest.py code)
for i in range(3):
    ax = plt.subplot(gs[1, i])
    ratio = np.loadtxt(f'dynamicalstate_senstitivitytest_b_pulse{i+1}_d90_ratio.txt')
    im = ax.imshow(ratio, cmap=cm, aspect='auto', origin='lower', extent=[0, (nb_time_flux-1)*dt, l_param[0], l_param[-1]])
    ax.set_xlabel('Time [d]', fontsize=10)
    ax.set_ylabel('b [mmolC m$^{-3}d^{-1}$]', fontsize=10)
    ax.axhline(y=0.08, color='gray', linestyle='--')
    ax.annotate(r'$P_2$', xy=(0.1, -0.49), xycoords='axes fraction', fontsize=12, ha='right')
    ax.annotate(r'$P_1$', xy=(0.9, -0.49), xycoords='axes fraction', fontsize=12, ha='left')
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.25, shrink=0.7, extend='both', ticks=[0.4, 0.5, 0.6])
    cbar.set_label(r'$R$ []', fontsize=10)
    im.set_clim(0.4, 0.6)
    ax.set_xlim(0, 90)
    
# Third line
for i in range(3):
    ax = plt.subplot(gs[2, i])
    ratio = np.loadtxt(f'dynamicalstate_senstitivitytest_b_pulse{i+1}_d90_Z.txt')
    im = ax.imshow(ratio, cmap=cm2, aspect='auto', origin='lower', extent=[0, (nb_time_flux-1)*dt, l_param[0], l_param[-1]])
    ax.set_xlabel('Time [d]', fontsize=10)
    ax.set_ylabel('b [mmolC m$^{-3}d^{-1}$]', fontsize=10)
    ax.axhline(y=0.08, color='gray', linestyle='--')
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.25, shrink=0.7)
    cbar.set_label(r'$Z$ [mmolC m$^{-3}$]', fontsize=10)
    im.set_clim(0.4, 0.6)
    ax.set_xlim(0, 90)
    
plt.tight_layout()
plt.show()
plt.savefig('../figures/dynamicaltest90d.pdf', format='pdf')