"""

Plot the dynaical-state outputs of the model as a function of a range of amplitude b. Create the figure 8 of the paper.

"""

import matplotlib.pyplot as plt
import numpy as npy
import numpy as np
import scienceplots
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
import sys
sys.path.append('../')

plt.style.use(['science','no-latex'])
plt.close('all')


## Configuration of the sensitivity test (identical to the configuration of the outpts_dynamicalstate_sensitivitytest.py code)
dt = 0.1
end_time_flux = 90
time_flux2 = np.arange(0, end_time_flux, dt)
nb_time = len(time_flux2)
min_param = 0.01
max_param = 0.1
l_param = np.linspace(min_param, max_param, 10)

## Color map configuration
colors = ['darkgreen', '#ffffff', 'chartreuse']
n_bins = 100
cmap_name = 'custom_green'
cm = plt.cm.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
plt.rc('font', size=7)


fig = plt.figure(figsize=(10, 7))
gs = GridSpec(3, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1, 1]) 

# Load data from the outputs_dynamicalstate.py code
data1 = np.loadtxt('dynamicalstate_pulse1_b0.08_d90.txt')
time_flux1, P1_1, P2_1, Z_1, PO4_1, Psupply1 = data1.T
data2 = np.loadtxt('dynamicalstate_pulse2_b0.08_d90.txt')
time_flux2, P1_2, P2_2, Z_2, PO4_2, Psupply2 = data2.T
data3 = np.loadtxt('dynamicalstate_pulse3_b0.08_d90.txt')
time_flux3, P1_3, P2_3, Z_3, PO4_3, Psupply3 = data3.T

# Loop to create the subplot

# First line
for i, data in enumerate([(time_flux1, P1_1, P2_1, Z_1, PO4_1, Psupply1),
                           (time_flux2, P1_2, P2_2, Z_2, PO4_2, Psupply2),
                           (time_flux3, P1_3, P2_3, Z_3, PO4_3, Psupply3)]):
    ax = plt.subplot(gs[0, i])
    time_flux, P1, P2, Z, PO4, Psupply = data
    ax.plot(time_flux, Psupply, label=r'$P_{supply}$', color="darkgray")
    ax.plot(time_flux, PO4, label=r'$PO4$', color="magenta")
    ax.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
    ax.plot(time_flux, P2, label=r'$P_2$', color="green")
    ax.plot(time_flux, Z, label=r'$Z$', color="aqua")
    ax.set_xlabel('Time [d]', fontsize=10)
    ax.set_ylabel('Masses [mmolCm$^{-3}$]', fontsize=10)
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
    im = ax.imshow(ratio, cmap=cm, aspect='auto', origin='lower', extent=[0, (nb_time-1)*dt, l_param[0], l_param[-1]])
    ax.set_xlabel('Time [d]', fontsize=10)
    ax.set_ylabel('b [mmolCm$^{-3}d^{-1}$]', fontsize=10)
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
    im = ax.imshow(ratio, cmap='Greys', aspect='auto', origin='lower', extent=[0, (nb_time-1)*dt, l_param[0], l_param[-1]])
    ax.set_xlabel('Time [d]', fontsize=10)
    ax.set_ylabel('b [mmolCm$^{-3}d^{-1}$]', fontsize=10)
    ax.axhline(y=0.08, color='gray', linestyle='--')
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.25, shrink=0.7)
    cbar.set_label(r'$Z$ [mmolCm$^{-3}$]', fontsize=10)
    im.set_clim(0.4, 0.6)
    ax.set_xlim(0, 90)
    
plt.tight_layout()
plt.show()
plt.savefig('../figures/dynamicaltest90d.pdf', format='pdf')