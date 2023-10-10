"""

Calculate and record the dynamical-state outputs of the model as a function of the value chosen for the amplitude b.

"""

from growth_model_2P1Z_v10 import growth_model_2P1Z_v10 
from fluxes import pulsedflux_stepfunction
import matplotlib.pyplot as plt
import numpy as npy
import numpy as np

# Define the set up
dt = 0.1
## To calculate steady-state equilibrium solutions
end_time = 2000
time = np.arange(1, end_time, dt)
Psupply_ini = 0.05
Psupply_cst =[Psupply_ini] * len(time)
## Pulse characteristics
b = 0.08
nbpulse = 3
## To calculate dynamical-state equilibrium solutions
end_time_flux = 90
time_flux = npy.arange(0,end_time_flux,dt)
nb_time=len(time_flux)


# Use the pulsedflux_stepfunction in the fluxes.py code in the "utils" directory
if nbpulse == 3 : 
    [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=b,C=Psupply_ini, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 25},{'t1': 35, 't2': 40}])
if nbpulse == 2 : 
    [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=b,C=Psupply_ini, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 25}])
if nbpulse == 1 : 
    [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=b,C=Psupply_ini, pulsations=[{'t1': 5, 't2': 10}])
    
    
# Steady-state solutions used to initiate dynamic state run
[P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply_cst,time,dt)
P1_new = P1[len(P1)-1]
P2_new = P2[len(P2)-1]
Z_new = Z[len(Z)-1]
PO4_new = PO4[len(PO4)-1]
ratio = P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1])

# Dynamical state solutions
[P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply,time_flux,dt,P1_ini=P1_new,P2_ini=P2_new,Z_ini=Z_new,PO4_ini=PO4_new)

Psupply_list = Psupply.tolist()
data_filename = f'../outputs/dynamicalstate_pulse{nbpulse}_b{b}_d{end_time_flux}.txt'
data = np.column_stack((time_flux, P1, P2, Z, PO4,Psupply_list))
np.savetxt(data_filename, data, header="time_flux P1 P2 Z PO4 Psupply")