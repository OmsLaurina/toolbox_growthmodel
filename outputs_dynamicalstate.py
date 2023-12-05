"""
* Calculate and record the dynamical-state outputs of the model as a function of the value chosen for the amplitude b.
* The configuration is defined within the "set_up" function in the "set_up.py" script.
"""

from growth_model import growth_model
from fluxes import pulsedflux_stepfunction
import numpy as np
from set_up import set_up_dynamicalstate, set_up_steadystate

# Configuration
dt, end_time, time, Psupply, Psupply_cst, Psupply_senstest, min_param, max_param,l_param = set_up_steadystate()
dt, Psupply_ini, b, nbpulse, end_time_flux, time_flux, nb_time_flux = set_up_dynamicalstate()

## To calculate dynamical-state equilibrium solutions

# Use the pulsedflux_stepfunction in the fluxes.py code in the "utils" directory
if nbpulse == 3 : 
    [Psupply_dyn,arg] = pulsedflux_stepfunction(time_flux, A=b,C=Psupply_ini, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 25},{'t1': 35, 't2': 40}])
if nbpulse == 2 : 
    [Psupply_dyn,arg] = pulsedflux_stepfunction(time_flux, A=b,C=Psupply_ini, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 25}])
if nbpulse == 1 : 
    [Psupply_dyn,arg] = pulsedflux_stepfunction(time_flux, A=b,C=Psupply_ini, pulsations=[{'t1': 5, 't2': 10}])  
    
# Steady-state solutions used to initiate dynamic state run
[P1,P2,Z,PO4,arg]=growth_model(Psupply,time,dt)
P1_new = P1[len(P1)-1]
P2_new = P2[len(P2)-1]
Z_new = Z[len(Z)-1]
PO4_new = PO4[len(PO4)-1]
ratio = P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1])

# Dynamical state solutions
[P1,P2,Z,PO4,arg]=growth_model(Psupply_dyn,time_flux,dt,P1_ini=P1_new,P2_ini=P2_new,Z_ini=Z_new,PO4_ini=PO4_new)

Psupply_list = Psupply_dyn.tolist()
data_filename = f'../outputs/dynamicalstate_pulse{nbpulse}_b{b}_d{end_time_flux}.txt'
data = np.column_stack((time_flux, P1, P2, Z, PO4,Psupply_list))
np.savetxt(data_filename, data, header="time_flux P1 P2 Z PO4 Psupply")