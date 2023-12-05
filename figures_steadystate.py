"""
* Plot the steady-state outputs of the model as a function of the value chosen for Psupply.
* The configuration is defined within the "set_up" function in the "set_up.py" script.
* First you need to create the outputs from the "outputs_steadystate.py" script.
"""

import matplotlib.pyplot as plt
import numpy as np
import scienceplots
from f_monod_hollingII import f_monod
from set_up import set_up_steadystate

plt.style.use(['science', 'no-latex'])
plt.close('all')

# Configuration
dt, end_time, time, Psupply, Psupply_cst, Psupply_senstest, min_param, max_param,l_param = set_up_steadystate()
test = Psupply_cst

# Choose the grazing control type
grazing = "diffgrazing"

# Load the data (create by "outputs_steadystate.py")
data = np.loadtxt(f'../outputs/steadystate_{grazing}_{test}.txt')
time, P1, P2, Z, PO4 = data.T

## Figure 1: Temporal evolution of biomass
plt.figure(1)
plt.rc('font', size=7)

plt.plot(time, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time, P2, label=r'$P_2$', color="green")
plt.plot(time, Z, label=r'$Z$', color="aqua")
plt.plot(time, PO4, label=r'$PO_4$', color="magenta")
plt.xlabel('Time [day]', fontsize=10)
plt.ylabel('Masses [mmolC.m\u207B\u00B3]', fontsize=10)
plt.xticks([0, 500, 1000, 1500, 2000])
legend = plt.legend(frameon=True, loc='upper right')
plt.gca().add_artist(legend)

figure_filename = f'../figures/steadystate_temporalevolution_{test}.pdf'
plt.savefig(figure_filename, format='pdf') 


## Figure 2:  Monod curves (use the f_monod function in the "utils" repertory)
# Obtain parameter values to calculate theoretical Monod curves
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

plt.figure(2)
plt.rc('font', size=7)

# Choose the theorical range of PO4
N_theo_PO4 = np.linspace(0, 6, len(time))

#Theorical Monod curves
f_monod(N_theo_PO4, arg['umax1'], arg['kP1'], "chartreuse")
f_monod(N_theo_PO4, arg['umax2'], arg['kP2'], "green")

#Modelled Monod curves
f_monod(PO4, arg['umax1'], arg['kP1'], "gray")
f_monod(PO4, arg['umax2'], arg['kP2'], "lightgray")

plt.xlabel(r'$PO_4$ [mmolC.m\u207B\u00B3]')
plt.ylabel(r'$\mu [d^{-1}]$')
plt.axvline(x=1.378, color='red', linestyle='--', label="[PO4] standard")
plt.axvline(x=PO4[-1], color='lightcoral', linestyle='--', label="[PO4] modelized")
plt.legend(handles=[plt.Line2D([0], [0], color='chartreuse', label=r'$P_1$'),
                    plt.Line2D([0], [0], color='green', label=r'$P_2$'),
                    plt.Line2D([0], [0], color='red', linestyle='--', label='[PO4] standard'),
                    plt.Line2D([0], [0], color='lightcoral', linestyle='--', label='[PO4] modelized')

                ])

figure_filename = f'../figures/Monodcurves_{test}.pdf'
plt.savefig(figure_filename, format='pdf')