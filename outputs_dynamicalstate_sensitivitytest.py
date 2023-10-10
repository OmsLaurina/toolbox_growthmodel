"""

Calculate and record the dynaical-state outputs of the model as a function of a range of amplitude b.

"""

import numpy as npy
import numpy as np
from growth_model_2P1Z_v10 import growth_model_2P1Z_v10 
from fluxes import pulsedflux_stepfunction


# Define the set up
dt = 0.1
## To calculate steady-state equilibrium solutions
end_time = 2000
time = np.arange(0, end_time, dt)
Psupply_ini = 0.05
Psupply_cst =[Psupply_ini] * len(time)
## Pulse characteristics
b = 0.08
nbpulse = 1
## To calculate dynamical-state equilibrium solutions
end_time_flux = 90
time_flux = npy.arange(0,end_time_flux,dt)
nb_time=len(time_flux)

# Steady-state solutions used to initiate dynamic state run
[P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply_cst, time, dt)
P1_eq = P1[len(P1)-1]
P2_eq = P2[len(P2)-1]
Z_eq = Z[len(Z)-1]
PO4_eq = PO4[len(PO4)-1]

# Number of value tested
n = 10

# Name of the parameters in the model you want to study 
param='b'

#Tested range of values
min_param = 0.01
max_param = 0.1
l_param = np.linspace(min_param,max_param,n)

# Create the empty matrices
P1_new = np.zeros((len(l_param), nb_time))
P2_new = np.zeros((len(l_param), nb_time))
PO4_new = np.zeros((len(l_param), nb_time))
Z_new = np.zeros((len(l_param), nb_time))
ratio = np.zeros((len(l_param),nb_time))
Psupply_new =  np.zeros((len(l_param), nb_time))
b_list=np.zeros(len(l_param))

filename=f"../outputs/dynamicalstate_senstitivitytest_{param}_pulse{nbpulse}_d{end_time_flux}.txt"

i=0

# amplitude range loop. Use the function pulsedflux_stepfunction in the fluxes.py code
for param in l_param:
    
    if nbpulse==1:
        [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=param,C=Psupply_ini, pulsations=[{'t1': 5, 't2': 10}])
        [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply,time_flux, dt, P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
    if nbpulse ==2:
           [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=param,C=Psupply_ini, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 25}])
           [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply,time_flux, dt, P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
    if nbpulse ==3:
            [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=param,C=Psupply_ini, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 25},{'t1': 35, 't2': 40}])
            [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply,time_flux, dt, P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
    
    P1_new[i, :] = np.array(P1)
    P2_new[i, :] = np.array(P2)
    PO4_new[i, :] = np.array(PO4)
    Z_new[i, :] = np.array(Z)
    ratio[i,:] = np.array(P1)/(np.array(P1)+ np.array(P2))
    Psupply_new[i, :] = np.array(Psupply)
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