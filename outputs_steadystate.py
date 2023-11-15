"""

Calculate and record the steady-state outputs of the model as a function of the value chosen for Psupply.

"""

import numpy as np
from growth_model import growth_model

# Define the set up
dt = 0.1
end_time = 2000
time = np.arange(0, end_time, dt)
Psupply_moy = 0.1
Psupply = [Psupply_moy] * len(time)

# For output file name (to save results obtained with different Psupply)
test = Psupply_moy

# Call the function
[P1, P2, Z, PO4, arg] = growth_model(Psupply, time, dt)

#### Testing the importance of differential grazing #####

grazing = 'diffgrazing'

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
prop_P1 = (P1[-1]/(P1[-1]+P2[-1]))*100
prop_P2 =  (P2[-1]/(P1[-1]+P2[-1]))*100

data_filename = f'../outputs/steadystate_{grazing}_{test}.txt'
data = np.column_stack((time, P1, P2, Z, PO4))
np.savetxt(data_filename, data, header="Time P1 P2 Z PO4")

