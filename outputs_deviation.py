"""

Calculate and record the steady-state outputs deviation from default values as a function of the value chosen for Psupply.

"""

import numpy as np
import matplotlib.pyplot as plt
from growth_model import growth_model
import scienceplots
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from set_up import set_up_steadystate
plt.style.use(['science','no-latex'])
plt.close('all')

dt, end_time, time, Psupply, Psupply_cst, Psupply_senstest, min_param, max_param,l_param = set_up_steadystate()
param_names = ['umax1', 'umax2', 'kP1', 'kP2', 'gmax1', 'gmax2', 'kZ1', 'kZ2']
param_names2 = [r'$u_{max,1}$', r'$u_{max,2}$', r'$K_{P,1}$', r'$K_{P,2}$', r'$g_{max,1}$', r'$g_{max,2}$', r'$K_{Z,1}$', r'$K_{Z,2}$']

# Default values of parameters
default_params = {
    'umax1': 1.9872,
    'umax2': 2.7648,
    'kP1': 1,
    'kP2': 3,
    'gmax1': 3.89,
    'gmax2': 0.43,
    'kZ1': 5,
    'kZ2': 20,
}

# Choose the facteur of deviation
tx = 2

#  Run 2: Default value /tx
run2 = [
    default_params['umax1']/tx,  # umax1
    default_params['umax2']/tx,  # umax2
    default_params['kP1']/tx,    # kP1
    default_params['kP2']/tx,    # kP2
    default_params['gmax1']/tx,  # gmax1
    default_params['gmax2']/tx,  # gmax2
    default_params['kZ1']/tx,    # kZ1
    default_params['kZ2']/tx,    # kZ2
]
params_run2 = {
    'umax1': run2[0],
    'umax2': run2[1],
    'kP1': run2[2],
    'kP2': run2[3],
    'gmax1': run2[4],
    'gmax2': run2[5],
    'kZ1': run2[6],
    'kZ2': run2[7]
}

# Run 3: Default value *tx
run3 = [
    default_params['umax1']*tx,  # umax1
    default_params['umax2']*tx,  # umax2
    default_params['kP1']*tx,    # kP1
    default_params['kP2']*tx,    # kP2
    default_params['gmax1']*tx,  # gmax1
    default_params['gmax2']*tx,  # gmax2
    default_params['kZ1']*tx,    # kZ1
    default_params['kZ2']*tx,    # kZ2
]
params_run3 = {
    'umax1': run3[0],
    'umax2': run3[1],
    'kP1': run3[2],
    'kP2': run3[3],
    'gmax1': run3[4],
    'gmax2': run3[5],
    'kZ1': run3[6],
    'kZ2': run3[7]
}


custom_colors =['#FFB6C1', '#FFD700', '#98FB98', '#87CEEB', '#FFA07A', '#9370DB', '#90EE90', '#F0E68C']
fig, axs = plt.subplots(5, 1, figsize=(10, 15))
variations = np.zeros((len(param_names),4))

# Define a threshold to force values close to 0 to be equal to 0
# Values greater than 0.002 are expected by the model. If they are lower, this means 0
threshold = 10**(-2)

legends = []

# Loop on parameter values
for i, param_name in enumerate(param_names):
    # Creating modified parameters for run2 and run3
    run2 = default_params.copy()
    run2[param_name] /= tx

    run3 = default_params.copy()
    run3[param_name] *= tx

    # Get solutions
    P1_default, P2_default, Z_default, PO4_default, _ = growth_model(Psupply, time, dt)
    P1_run2, P2_run2, Z_run2, PO4_run2, _ = growth_model(Psupply, time, dt, **run2)
    P1_run3, P2_run3, Z_run3, PO4_run3, _ = growth_model(Psupply, time, dt, **run3)
    R_default = np.array(P1_default) / (np.array(P1_default) + np.array(P2_default))
    R_run3 = np.array(P1_run3) / (np.array(P1_run3) + np.array(P2_run3))
    R_run2 = np.array(P1_run2) / (np.array(P1_run2) + np.array(P2_run2))
    
    # Graphically check solutions for each state variable and parameter value
    color = custom_colors[i % len(custom_colors)]
    
    axs[0].plot(time, PO4_run2, color=color)
    axs[0].plot(time, PO4_run3, linestyle='dotted',color=color)
    axs[0].axhline(y=PO4_default[-1], color='red', linestyle='--', label='Default')
    axs[0].set_title(r'$PO_4$')
    
    axs[1].plot(time, P1_run2, color=color)
    axs[1].plot(time, P1_run3, linestyle='dotted', color=color)
    axs[1].axhline(y=P1_default[-1], color='red', linestyle='--', label='Default')
    axs[1].set_title(r'$P1$')
    
    axs[2].plot(time, P2_run2, color=color)
    axs[2].plot(time, P2_run3,  linestyle='dotted',color=color)
    axs[2].axhline(y=P2_default[-1], color='red', linestyle='--', label='Default')
    axs[2].set_title(r'$P2$')
    
    axs[3].plot(time, Z_run2, label='Run2', color=color)
    axs[3].plot(time, Z_run3, label='Run3', linestyle='dotted',color=color)
    axs[3].axhline(y=Z_default[-1], color='red', linestyle='--', label='Default')
    axs[3].set_title('Z')
    
    axs[4].plot(time, R_run2, color=color)
    axs[4].plot(time, R_run3, linestyle='dotted',color=color)
    axs[4].axhline(y=R_default[-1], color='red', linestyle='--', label='Default')
    axs[4].set_title('R')
    
    # Get the default
    P1_default = np.mean(P1_default[-200:])
    P2_default = np.mean(P2_default[-200:])
    Z_default = np.mean(Z_default[-200:])
    PO4_default = np.mean(PO4_default[-200:])
    R_default = np.mean(R_default[-200:])
    if P2_default<threshold:
        P2_default = 0
        R_default = 1
    
    # Get the run2
    P1_run2 = np.mean(P1_run2[-200:])
    P2_run2 = np.mean(P2_run2[-200:])
    Z_run2 = np.mean(Z_run2[-200:])
    PO4_run2 = np.mean(PO4_run2[-200:])
    R_run2 = np.mean(R_run2[-200:])
    if P2_run2<threshold:
        P2_run2 = 0
        R_run2 = 1
    
    # Get the run3
    P1_run3 = np.mean(P1_run3[-200:])
    P2_run3 = np.mean(P2_run3[-200:])
    Z_run3 = np.mean(Z_run3[-200:])
    PO4_run3 = np.mean(PO4_run3[-200:])
    R_run3 = np.mean(R_run3[-200:])
    if P2_run3<threshold:
        P2_run3 = 0
        R_run3 = 1
    
    # Calculation of the variation
    P1_variation = ((P1_run3 - P1_run2) / P1_default) * 100
    if P2_default == 0:
        P2_variation = float('NaN')
    else:
        P2_variation = ((P2_run3 - P2_run2) / P2_default) * 100
    Z_variation = ((Z_run3 - Z_run2) / Z_default) * 100
    PO4_variation = ((PO4_run3 - PO4_run2) / PO4_default) * 100
    R_variation = ((R_run3 - R_run2) / R_default) * 100
    print(R_variation)
    
    variations[i,0:]= np.array([PO4_variation,P1_variation,P2_variation,R_variation])

    
plt.subplots_adjust(hspace=0.7)

for i, ax in enumerate(axs):
    ax.set_xlabel('Time [days]', fontsize = 9)
    ax.set_ylabel(r'Masse $[mmolC m^{-3}]$', fontsize = 9)

legend_patches = [mpatches.Patch(color=color, label=param_names2) for color, param_names2 in zip(custom_colors, param_names2)]
fig.legend(handles=legend_patches, loc='center right', bbox_to_anchor=(1.02, 0.5), frameon=False)

custom_lines = [Line2D([0], [0], color='black', linestyle='dotted', label='twice (run 3)'),
                Line2D([0], [0], color='black', linestyle='-', label='half (run 2)')]
fig.legend(handles=custom_lines, loc='center right', bbox_to_anchor=(1.02, 0.38), frameon=False)

plt.show()
plt.savefig(f'../figures/sensitivity_test_{Psupply_cst}.pdf', format='pdf')
np.savetxt(f'../outputs/sensitivity_variations_{Psupply_cst}.txt', variations, fmt='%.8f')