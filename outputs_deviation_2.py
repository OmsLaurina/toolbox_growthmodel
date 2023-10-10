#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 07:02:39 2023

@author: loms
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 15:32:50 2023

@author: loms
"""
import numpy as np
import numpy as npy
import matplotlib.pyplot as plt
from growth_model_2P1Z_v10 import growth_model_2P1Z_v10
import scienceplots
import pandas as pd
plt.style.use(['science','no-latex'])
plt.close('all')

dt = 0.1
end_time = 2000
time = npy.arange(0, end_time, dt)
p = 0.01
test=p
Psupply = [p] * len(time)

param_names = ['umax1', 'umax2', 'kP1', 'kP2', 'gmax1', 'gmax2', 'kZ1', 'kZ2']

# Default values
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

tx = 2

# Default value /tx
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

# Default value *tx
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

fig, axs = plt.subplots(5, 1, figsize=(10, 15))
custom_colors =['#FFB6C1', '#FFD700', '#98FB98', '#87CEEB', '#FFA07A', '#9370DB', '#90EE90', '#F0E68C']

variations = np.zeros((len(param_names),5))
threshold = 10**(-2)

# Boucle pour chaque paramètre
for i, param_name in enumerate(param_names):
    # Création de paramètres modifiés pour run2 et run3
    run2 = default_params.copy()
    run2[param_name] /= tx

    run3 = default_params.copy()
    run3[param_name] *= tx

    # Get solutions
    P1_default, P2_default, Z_default, PO4_default, _ = growth_model_2P1Z_v10(Psupply, time, dt)
    P1_run2, P2_run2, Z_run2, PO4_run2, _ = growth_model_2P1Z_v10(Psupply, time, dt, **run2)
    P1_run3, P2_run3, Z_run3, PO4_run3, _ = growth_model_2P1Z_v10(Psupply, time, dt, **run3)
    R_default = np.array(P1_default) / (np.array(P1_default) + np.array(P2_default))
    R_run3 = np.array(P1_run3) / (np.array(P1_run3) + np.array(P2_run3))
    R_run2 = np.array(P1_run2) / (np.array(P1_run2) + np.array(P2_run2))
    
    color = custom_colors[i % len(custom_colors)]
    axs[0].plot(time, P1_run2, color=color)
    axs[0].plot(time, P1_run3, linestyle='dotted', color=color)
    axs[0].axhline(y=P1_default[-1], color='red', linestyle='--', label='Default')
    axs[0].set_title('P1')
    # axs[0].legend()
    
    axs[1].plot(time, P2_run2, color=color)
    axs[1].plot(time, P2_run3,  linestyle='dotted',color=color)
    axs[1].axhline(y=P2_default[-1], color='red', linestyle='--', label='Default')
    axs[1].set_title('P2')
    # axs[1].legend()
    
    axs[2].plot(time, Z_run2, label='Run2', color=color)
    axs[2].plot(time, Z_run3, label='Run3', linestyle='dotted',color=color)
    axs[2].axhline(y=Z_default[-1], color='red', linestyle='--', label='Default')
    axs[2].set_title('Z')
    # axs[2].legend()
    
    axs[3].plot(time, PO4_run2, color=color)
    axs[3].plot(time, PO4_run3, linestyle='dotted',color=color)
    axs[3].axhline(y=PO4_default[-1], color='red', linestyle='--', label='Default')
    axs[3].set_title('PO4')
    # axs[3].legend()
    
    axs[4].plot(time, R_run2, color=color)
    axs[4].plot(time, R_run3, linestyle='dotted',color=color)
    axs[4].axhline(y=R_default[-1], color='red', linestyle='--', label='Default')
    axs[4].set_title('R')
    # axs[4].legend()
    
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
    
    variations[i,0:]= np.array([P1_variation,P2_variation,Z_variation,PO4_variation,R_variation])

plt.subplots_adjust(hspace=0.5)
plt.show()
plt.savefig(f'../outputs/sensitivity_test_{p}.pdf', format='pdf')

# Save data to a text file
np.savetxt(f'../outputs/sensitivity_variations_{test}.txt', variations, fmt='%.8f')