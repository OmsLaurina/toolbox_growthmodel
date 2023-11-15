"""

Calculate and record the coexistence threshold deviation from default values as a function of the value chosen for Psupply.

"""

import numpy as np
import numpy as npy
import matplotlib.pyplot as plt
from growth_model import growth_model
import scienceplots
import pandas as pd
plt.style.use(['science','no-latex'])
plt.close('all')

# Define the set up
dt = 0.1
end_time = 2000
time = npy.arange(0, end_time, dt)
param_names = ['umax1', 'umax2', 'kP1', 'kP2', 'gmax1', 'gmax2', 'kZ1', 'kZ2']

# To get the coexistence threshold 
Psupply_moy = 1
Psupply_list = [Psupply_moy]*len(time)
Psupply_arr = np.array(Psupply_list)
param = 'Psupp'
# Tested range of values
min_param = 0.01
max_param = 0.096
l_param = np.linspace(min_param,max_param,100)

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

# Préallouer des tableaux pour stocker les résultats
R_default_2 = np.zeros(len(l_param))
R_run2_2 = np.zeros(len(l_param))
R_run3_2 = np.zeros(len(l_param))

# # Listes pour stocker les valeurs de param
# coex_threshold_default_values = []
# coex_threshold_run2_values = []
# coex_threshold_run3_values = []
# coex_threshold_variations = []

# Variable de contrôle pour suivre combien de seuils ont été trouvés
# thresholds_found = 0

fig, axs = plt.subplots(2, 4, figsize=(12, 6))
axs = axs.flatten()
                        
# Loop on parameter values
for i, param_name in enumerate(param_names):
    
    # Creating modified parameters for run2 and run3
    run2 = default_params.copy()
    run2[param_name] /= tx

    run3 = default_params.copy()
    run3[param_name] *= tx

    # Variables pour stocker la valeur de param pour chaque run
    coex_threshold_default = np.nan
    coex_threshold_run2 = np.nan
    coex_threshold_run3 = np.nan
    
    bo = True

    # Get solutions for the current param value
    for j, param in enumerate(l_param):
        #print(param)
        Psupply_value = Psupply_arr * param
        P1_default, P2_default, Z_default, PO4_default, _ = growth_model(Psupply_value, time, dt)
        P1_run2, P2_run2, Z_run2, PO4_run2, _ = growth_model(Psupply_value, time, dt, **run2)
        P1_run3, P2_run3, Z_run3, PO4_run3, _ = growth_model(Psupply_value, time, dt, **run3)
        R_default = np.array(P1_default) / (np.array(P1_default) + np.array(P2_default))
        R_run2 = np.array(P1_run2) / (np.array(P1_run2) + np.array(P2_run2))
        R_run3 = np.array(P1_run3) / (np.array(P1_run3) + np.array(P2_run3))

        # Calculate the means of R values for the last 200 time steps
        R_default_2[j] = np.mean(R_default[-2:])
        R_run2_2[j] = np.mean(R_run2[-2:])
        R_run3_2[j] = np.mean(R_run3[-2:])
        
    #     # Check if R becomes less than 1 for the first time and store the corresponding param value
    #     if R_default_2[j] < 1 and bo:
    #         coex_threshold_default = param
    #         coex_threshold_default_values.append(param)
    #         thresholds_found += 1

    #     if R_run2_2[j] < 1 and bo:
    #         coex_threshold_run2 = param
    #         coex_threshold_run2_values.append(param)
    #         thresholds_found += 1

    #     if R_run3_2[j] < 1 and bo:
    #         coex_threshold_run3 = param
    #         coex_threshold_run3_values.append(param)
    #         thresholds_found += 1
        
    #     # Check if the 3 threshold are found
    #     if thresholds_found == 3 and bo:
    #         bo = False
    
    # # Calcul du pourcentage de variation du seuil de coexistence pour l'itération actuelle
    # coex_threshold_variation = ((coex_threshold_run3 - coex_threshold_run2) / coex_threshold_default) * 100
    # coex_threshold_variations.append(coex_threshold_variation)
    # np.savetxt('../outputs/sensitivity_coexthreshold.txt', coex_threshold_variations, fmt='%.8f')
    
    print(l_param[9])
    
    # Plot R values for the current param_name in the i-th subplot
    axs[i].plot(l_param, R_default_2, label='R_default', color='blue')
    axs[i].plot(l_param, R_run2_2, label='R_run2', color='green', linestyle='--')
    axs[i].plot(l_param, R_run3_2, label='R_run3', color='red', linestyle='-.')

    axs[i].set_ylabel('R')
    axs[i].set_title(f'R vs. {param_name}')
    axs[i].grid(True)

axs[-1].set_xlabel(param_name)
axs[-1].legend()

plt.tight_layout()
plt.show()
plt.savefig("../figures/deviation_coexthreshold.pdf", format='pdf')