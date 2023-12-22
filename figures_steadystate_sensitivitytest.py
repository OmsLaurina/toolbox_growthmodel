"""
* Plot the steady-state outputs of the model as a function of a range of Psupply value. Plot the figure 5 and 6 of the paper.
* The configuration is defined within the "set_up" function in the "set_up.py" script.
* First you need to create the outputs from the "outputs_steadystate_sensitivitytest.py" script.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scienceplots
from f_monod_hollingII import f_monod
from set_up import set_up_steadystate
from scipy.signal import savgol_filter

plt.style.use(['science','no-latex'])
plt.close('all')

# Configuration
dt, end_time, time, Psupply, Psupply_cst, Psupply_senstest, min_param, max_param,l_param = set_up_steadystate()

# Choose the grazing control type
grazing = 'diffgrazing'

# axis name
name_param = r'$P_{\mathrm{supply}}$ [mmolC m$^{-3}d^{-1}$]'

### Figure 1: Effect of Psupply on the R ratio and time trend for the two extreme Psupply values
name_pdf = f"../figures/{grazing}.pdf"

fig = plt.figure(figsize=(12,3))
gs = GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.3)

# Load data from outputs_steadystate_sensitivitytest.py code
P_1 = np.loadtxt(f'steadystate_senstitivitytest_Psupp_{grazing}_P1.txt')
P_2 = np.loadtxt(f'steadystate_senstitivitytest_Psupp_{grazing}_P2.txt')
ratio = np.loadtxt(f'steadystate_senstitivitytest_Psupp_{grazing}_ratio.txt')
P_O4 = np.loadtxt(f'steadystate_senstitivitytest_Psupp_{grazing}_PO4.txt')
Z_ = np.loadtxt(f'steadystate_senstitivitytest_Psupp_{grazing}_Z.txt')

# Middle panel
ax = fig.add_subplot(gs[1])
sc = ax.scatter(l_param, ratio, c=P_1+P_2, cmap='Wistia')
# sc = ax.scatter(l_param, Z_, c=P_1 + P_2, cmap='Wistia')
ax.plot(l_param, ratio, color='black')
ax.axvline(x=0.01, color='gray', linestyle='--',linewidth=1.5)
ax.axvline(x=0.02, color='gray', linestyle=':',linewidth=1.5)
ax.axvline(x=0.1, color='gray', linestyle='--',linewidth=1.5)
if grazing == 'diffgrazing':
    ax.axvline(x=0.05, color='red', linestyle='--')
ax.set_xlabel(name_param, fontsize=10)
ax.set_ylabel('R []', fontsize=10)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
divider = make_axes_locatable(ax)
cax = divider.append_axes("top", size="5%", pad=0.2)  # Ajustez la taille et le pad selon vos besoins
cbar = fig.colorbar(sc, cax=cax, orientation="horizontal", ticks=MaxNLocator(nbins=3))
cbar.ax.xaxis.set_label_position('top')
cbar.set_label(r'$P_1+P_2$ [mmolC m$^{-3}$]', labelpad=2)

# Left panel
# Load data from outputs_steadystate.py code
test = 0.01 
data = np.loadtxt(f'../outputs/steadystate_{grazing}_{test}.txt')
time, P1, P2, Z, PO4 = data.T
ax = fig.add_subplot(gs[0])
ax.plot(time, PO4, label=r'$PO_4$', color="magenta")
ax.plot(time, P1, label=r'$P_1$', color="chartreuse")
ax.plot(time, P2, label=r'$P_2$', color="green")
ax.plot(time, Z, label=r'$Z$', color="aqua")
ax.set_xlabel('Time [d]', fontsize=10)
ax.set_ylabel(r'Masses [mmolC m$^{-3}$]', fontsize=10)
ax.set_xticks([0, 500, 1000])
ax.set_xlim(0,1000)
legend = ax.legend(frameon=False, loc='upper right')
for label in legend.get_texts():
    label.set_fontsize(7)

# Right panel    
# Load data from outputs_steadystate.py code
test = 0.1
data = np.loadtxt(f'../outputs/steadystate_{grazing}_{test}.txt')
time, P1, P2, Z, PO4 = data.T
ax = fig.add_subplot(gs[2])
ax.plot(time, PO4, label=r'$PO_4$', color="magenta")
ax.plot(time, P1, label=r'$P_1$', color="chartreuse")
ax.plot(time, P2, label=r'$P_2$', color="green")
ax.plot(time, Z, label=r'$Z$', color="aqua")
ax.set_xlabel('Time [d]', fontsize=10)
ax.set_ylabel(r'Masses [mmolC m$^{-3}$]', fontsize=10)
ax.set_xticks([0, 500, 1000])
ax.set_xlim(0,1000)
    
plt.tight_layout()
plt.savefig(name_pdf, format='pdf')
plt.show()

### Figure 2: Only the R-ratio
fig = plt.figure(2)

sc = plt.scatter(l_param, ratio, c=P_1+P_2, cmap='Wistia')
plt.plot(l_param, ratio, color='black')
# plt.axvline(x=0.01, color='gray', linestyle='--', linewidth=1.5)
plt.axvline(x=0.03, color='gray', linestyle=':', linewidth=1.5)
# plt.axvline(x=0.096, color='gray', linestyle='--', linewidth=1.5)

if grazing == 'diffgrazing':
    plt.axvline(x=0.05, color='red', linestyle='--')
    # plt.axhline(y=0.5, color='red', linestyle='--')

plt.xlabel(name_param, fontsize=10)
plt.ylabel('R []', fontsize=10)
cbar = plt.colorbar(sc, ticks=MaxNLocator(nbins=3))
cbar.ax.xaxis.set_label_position('top')
cbar.set_label(r'$P_1+P_2$ [mmolC m$^{-3}$]', labelpad=2)
plt.yticks([0, 0.25, 0.5, 0.75, 1])

index_ratio_05 = np.argmax(ratio <= 0.5)
l_param_at_ratio_05 = l_param[index_ratio_05]
print("Index where ratio is lower than 0.5", index_ratio_05)
print("l_param corresponding value :", l_param_at_ratio_05)

# # Determine the inflexion point
# smoothed_curve = savgol_filter(ratio, window_length=5, polyorder=2)
# smoothed_curve = smoothed_curve[40:70]
# l_param = l_param[40:70]
# inflexion_point_index = np.argmax(np.gradient(np.gradient(smoothed_curve)) > 0)
# inflexion_point = (l_param[inflexion_point_index], smoothed_curve[inflexion_point_index])
# plt.annotate(f'Inflexion Point\n{name_param}={inflexion_point[0]:.2f}, R={inflexion_point[1]:.2f}',
#              xy=inflexion_point, xytext=(0.05, 0.75),
#              arrowprops=dict(facecolor='black', shrink=0.05))


plt.savefig(f'../figures/R-ratio_{grazing}.pdf', format='pdf')


# Figure 3: Monod curve with the range of PO4 cover by the model
plt.figure(3)
plt.rc('font', size=7)
name_pdf_monod = '../figures/monod.pdf'
arg = {
    'umax1': 1.9872,
    'umax2': 2.7648,
    'kP1': 1,
    'kP2': 3,
    'gmax1':3.89,
    'gmax2':0.43,
    'kZ1': 5,
    'kZ2': 20,
}
N_theo_PO4 = np.linspace(0, 6, len(time))
f_monod(N_theo_PO4/130, arg['umax1'], arg['kP1']/130,"chartreuse")
f_monod(N_theo_PO4/130, arg['umax2'], arg['kP2']/130,"green")
# f_monod(P_O4, arg['umax1'], arg['kP1'],"gray")
# f_monod(P_O4, arg['umax2'], arg['kP2'],"lightgray")
plt.xlabel(r'$PO_4 [mmolP\,m^{-3}]$')
plt.ylabel(r'$\mu\ [d^{-1}]$')
plt.ylim(0,2)
plt.axvline(x=0.013, color='red', linestyle='--', label="[PO4] observed")
plt.axvline(x=min(P_O4)/130, color='lightcoral', linestyle='--', label="Min PO4")
plt.axvline(x=max(P_O4)/130, color='lightcoral', linestyle='--', label="Max PO4")
plt.fill_betweenx(np.linspace(0, max(plt.ylim()), len(l_param)), min(P_O4)/130, max(P_O4)/130, color='red', alpha=0.3, label='Range')

plt.legend(handles=[plt.Line2D([0], [0], color='chartreuse', label=r'$P_1$'),
                plt.Line2D([0], [0], color='green', label=r'$P_2$'),
                plt.Line2D([0], [0], color='red', linestyle='--', label=r'Observed [$PO_4$]'),
                plt.Line2D([0], [0], color='lightcoral', linestyle='--', label=r'Modeled [$PO_4$] range')
                ])

plt.savefig(name_pdf_monod, format='pdf')