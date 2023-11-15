"""

plot the steady-state outputs of the model as a function of a range of Psupply value. Plot the figure 5 and 6 of the paper.

"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scienceplots
from f_monod_hollingII import f_monod

plt.style.use(['science','no-latex'])
plt.close('all')

## Configuration of the Psupply sensitivity test (identical to the configuration of the outpts_steadystate_sensitivitytest.py code)
name_param = r'$P_{\mathrm{supply}}$ [mmolCm$^{-3}d^{-1}$]'
n = 100
min_param = 0.01
max_param = 0.1
l_param = np.linspace(min_param, max_param, n)
grazing = "diffgrazing"

# Figure 1 (figure 6 in article): Effect of Psupply on the R ratio and time trend for the two extreme Psupply values.
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
cbar.set_label(r'$P_1+P_2$ [mmolCm$^{-3}$]', labelpad=2)


# ratio_min = 0.25
# ratio_max = 0.3
# l_param_within_range = []
# for i in range(len(l_param)):
#     current_ratio = ratio[i]
#     if ratio_min <= current_ratio <= ratio_max:
#         l_param_within_range.append(l_param[i])

# Left panel
# Load data from outputs_steadystate.py code
test = 0.01 
data = np.loadtxt(f'../outputs/steadystate_{grazing}_{test}.txt')
time, P1, P2, Z, PO4 = data.T
ax = fig.add_subplot(gs[0])
ax.plot(time, P1, label=r'$P_1$', color="chartreuse")
ax.plot(time, P2, label=r'$P_2$', color="green")
ax.plot(time, Z, label=r'$Z$', color="aqua")
ax.plot(time, PO4, label=r'$PO4$', color="magenta")
ax.set_xlabel('Time [d]', fontsize=10)
ax.set_ylabel(r'Masses [mmolCm$^{-3}$]', fontsize=10)
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
ax.plot(time, P1, label=r'$P_1$', color="chartreuse")
ax.plot(time, P2, label=r'$P_2$', color="green")
ax.plot(time, Z, label=r'$Z$', color="aqua")
ax.plot(time, PO4, label=r'$PO4$', color="magenta")
ax.set_xlabel('Time [d]', fontsize=10)
ax.set_ylabel(r'Masses [mmolCm$^{-3}$]', fontsize=10)
ax.set_xticks([0, 500, 1000])
ax.set_xlim(0,1000)
    
plt.tight_layout()
plt.savefig(name_pdf, format='pdf')
plt.show()

# Only the R-ratio
fig = plt.figure(2)

sc = plt.scatter(l_param, ratio, c=P_1+P_2, cmap='Wistia')
plt.plot(l_param, ratio, color='black')
# plt.axvline(x=0.01, color='gray', linestyle='--', linewidth=1.5)
# plt.axvline(x=0.02, color='gray', linestyle=':', linewidth=1.5)
# plt.axvline(x=0.096, color='gray', linestyle='--', linewidth=1.5)

if grazing == 'diffgrazing':
    plt.axvline(x=0.05, color='red', linestyle='--')

plt.xlabel(name_param, fontsize=10)
plt.ylabel('R []', fontsize=10)
cbar = plt.colorbar(sc, ticks=MaxNLocator(nbins=3))
cbar.ax.xaxis.set_label_position('top')
cbar.set_label(r'$P_1+P_2$ [mmolCm$^{-3}$]', labelpad=2)
plt.yticks([0, 0.25, 0.5, 0.75, 1])
plt.savefig(f'../figures/R-ratio_{grazing}.pdf', format='pdf')


# # Figure 2 (figure 5 in article): Effect of Psupply on the R ratio and time trend for the two extreme Psupply values.
# plt.figure(3)
# plt.rc('font', size=7)
# name_pdf_monod = '../figures/monod.pdf'
# arg = {
#     'umax1': 1.9872,
#     'umax2': 2.7648,
#     'kP1': 1,
#     'kP2': 3,
#     'gmax1':3.89,
#     'gmax2':0.43,
#     'kZ1': 5,
#     'kZ2': 20,
# }
# N_theo_PO4 = np.linspace(0, 6, len(time))
# f_monod(N_theo_PO4, arg['umax1'], arg['kP1'],"chartreuse")
# f_monod(N_theo_PO4, arg['umax2'], arg['kP2'],"green")
# # f_monod(P_O4, arg['umax1'], arg['kP1'],"gray")
# # f_monod(P_O4, arg['umax2'], arg['kP2'],"lightgray")
# plt.xlabel('PO4 [mmolCm\u207B\u00B3]')
# plt.ylabel(r'$\mu\ [d^{-1}]$')
# plt.ylim(0,2)
# plt.axvline(x=1.378, color='red', linestyle='--', label="[PO4] observed")
# plt.axvline(x=min(P_O4), color='lightcoral', linestyle='--', label="Min PO4")
# plt.axvline(x=max(P_O4), color='lightcoral', linestyle='--', label="Max PO4")
# plt.fill_betweenx(np.linspace(0, max(plt.ylim()), n), min(P_O4), max(P_O4), color='red', alpha=0.3, label='Range')

# plt.legend(handles=[plt.Line2D([0], [0], color='chartreuse', label=r'$P_1$'),
#                 plt.Line2D([0], [0], color='green', label=r'$P_2$'),
#                 plt.Line2D([0], [0], color='red', linestyle='--', label='Observed [PO4]'),
#                 plt.Line2D([0], [0], color='lightcoral', linestyle='--', label='Modeled [PO4] range')
#                 ])

# plt.savefig(name_pdf_monod, format='pdf')