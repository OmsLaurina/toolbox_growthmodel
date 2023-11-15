"""

Calculate and record the dynamical-state outputs of the model as a function of a range of amplitude b.

"""

import numpy as np
from growth_model import growth_model
from fluxes import pulsedflux_stepfunction
from set_up import set_up_steadystate

# Define the set up

dt, end_time, time, Psupply, Psupply_cst, Psupply_senstest, l_param = set_up_steadystate() # modify it in the set_up.py script

## Pulse characteristics
nbpulse = 1

## To calculate dynamical-state equilibrium solutions

end_time_flux = 90
time_flux = np.arange(0,end_time_flux,dt)
nb_time=len(time_flux)

# Steady-state solutions used to initiate dynamic state run
[P1,P2,Z,PO4,arg]=growth_model(Psupply, time, dt)
P1_eq = P1[len(P1)-1]
P2_eq = P2[len(P2)-1]
Z_eq = Z[len(Z)-1]
PO4_eq = PO4[len(PO4)-1]

# Number of value tested
n = 100

# Name of the parameters in the model you want to study 
param= 'b'

#Tested range of values
min_param = 0.01
max_param = 0.1
l_param2 = np.linspace(min_param,max_param,n)

# Create the empty matrices
P1_new = np.zeros((len(l_param), nb_time))
P2_new = np.zeros((len(l_param), nb_time))
PO4_new = np.zeros((len(l_param), nb_time))
Z_new = np.zeros((len(l_param), nb_time))
ratio = np.zeros((len(l_param),nb_time))
Psupply_new =  np.zeros((len(l_param), nb_time))

filename=f"../outputs/dynamicalstate_senstitivitytest_{param}_pulse{nbpulse}_d{end_time_flux}.txt"

i=0

# amplitude range loop. Use the function pulsedflux_stepfunction in the fluxes.py code
for param in l_param2:
    
    if nbpulse==1:
        [Psupply_dyn,arg] = pulsedflux_stepfunction(time_flux, A=param,C=Psupply_cst, pulsations=[{'t1': 5, 't2': 10}])
        [P1,P2,Z,PO4,arg]=growth_model(Psupply_dyn,time_flux, dt, P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
    if nbpulse ==2:
           [Psupply_dyn,arg] = pulsedflux_stepfunction(time_flux, A=param,C=Psupply_cst, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 25}])
           [P1,P2,Z,PO4,arg]=growth_model(Psupply_dyn,time_flux, dt, P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
    if nbpulse ==3:
            [Psupply_dyn,arg] = pulsedflux_stepfunction(time_flux, A=param,C=Psupply_cst, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 25},{'t1': 35, 't2': 40}])
            [P1,P2,Z,PO4,arg]=growth_model(Psupply_dyn,time_flux, dt, P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
    
    P1_new[i, :] = np.array(P1)
    P2_new[i, :] = np.array(P2)
    PO4_new[i, :] = np.array(PO4)
    Z_new[i, :] = np.array(Z)
    ratio[i,:] = np.array(P1)/(np.array(P1)+ np.array(P2))
    Psupply_new[i, :] = np.array(Psupply_dyn)
    i=i+1
    print(i)
       
filename_P1 = filename.replace('.txt', '_P1.txt')
filename_P2 = filename.replace('.txt', '_P2.txt')
filename_ratio = filename.replace('.txt', '_ratio.txt')
filename_PO4 = filename.replace('.txt', '_PO4.txt')
filename_Z = filename.replace('.txt', '_Z.txt')
filename_Psupply = filename.replace('.txt', '_Psupply.txt')
        
np.savetxt(filename_P1, P1_new)
np.savetxt(filename_P2, P2_new)
np.savetxt(filename_ratio, ratio)
np.savetxt(filename_PO4, PO4_new)
np.savetxt(filename_Z, Z_new)
np.savetxt(filename_Psupply, Psupply_new)