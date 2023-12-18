"""
* Calculate and record the steady-state outputs of the model as a function of the value chosen for Psupply.
* The configuration is defined within the "set_up" function in the "set_up.py" script.
"""

import numpy as np
from growth_model import growth_model
from set_up import set_up_steadystate
import matplotlib.pyplot as plt

# Configuration
dt, end_time, time, Psupply, Psupply_cst, Psupply_senstest, min_param, max_param,l_param = set_up_steadystate()

# For output file name (to save results obtained with different Psupply)
test = Psupply_cst

# Call the growth function
[P1, P2, Z, PO4, arg] = growth_model(Psupply, time, dt)

#### Testing the importance of differential grazing controls #####

grazing = 'nograzing'

[P1, P2, Z, PO4, arg] = growth_model(Psupply, time, dt, gmax1=0, gmax2=0)
data_filename = f'../outputs/steadystate_{grazing}_{test}.txt'
data = np.column_stack((time, P1, P2, Z, PO4))
np.savetxt(data_filename, data, header="Time P1 P2 Z PO4")

###

grazing = 'equalgrazing'

[P1, P2, Z, PO4, arg] = growth_model(Psupply, time, dt, kZ2=5, gmax2=3.89)
data_filename = f'../outputs/steadystate_{grazing}_{test}.txt'
data = np.column_stack((time, P1, P2, Z, PO4))
np.savetxt(data_filename, data, header="Time P1 P2 Z PO4")

###

grazing = 'diffgrazing'

[P1, P2, Z, PO4, arg] = growth_model(Psupply, time, dt)

dP1 = []
dP2 = []
dPO4 = []
dZ = []

nb_time=len(time)
for i in range(0,nb_time-1):
    dP1next = (P1[i+1]-P1[i])/dt
    dP2next = (P2[i+1]-P2[i])/dt
    dPO4next = (PO4[i+1]-PO4[i])/dt
    dZnext = (Z[i+1]-Z[i])/dt
    
    dP1.append(dP1next)
    dP2.append(dP2next)
    dPO4.append(dPO4next)
    dZ.append(dZnext)
    
plt.plot(time[0:len(time)-1],dP1)
    
    
data_filename = f'../outputs/steadystate_{grazing}_{test}.txt'
data = np.column_stack((time, P1, P2, Z, PO4))
np.savetxt(data_filename, data, header="Time P1 P2 Z PO4")