"""

Calculate and record the steady-state outputs of the model as a function of a range of Psupply value.

"""
import numpy as npy
import numpy as np
from growth_model_2P1Z_v10 import growth_model_2P1Z_v10
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science', 'no-latex'])
plt.close('all')

# Define the set up
dt = 0.1
end_time = 2000
time = npy.arange(0,end_time,dt)
nb_time = len(time)
grazing = 'diffgrazing'

# Set Psupply equal to 1 to test its value
Psupply_moy = 1
Psupply = [Psupply_moy]*len(time)
Psupply_arr = np.array(Psupply)

# Number of tested values
n = 100

# Name of the parameters in the model you want to study 
param = 'Psupp' #for file name
name_param =r'$P_{supply}$' #for plot

# Tested range of values
min_param = 0.01
max_param = 0.096
l_param = np.linspace(min_param,max_param,n)
    
filename = f"../outputs/steadystate_senstitivitytest_{param}_{grazing}.txt"

# Create the empty matrices
P_1=np.zeros((len(l_param)))
P_2=np.zeros((len(l_param)))
ratio=np.zeros((len(l_param)))
P_O4=np.zeros((len(l_param)))
Z_=np.zeros((len(l_param)))
prop_P1=np.zeros((len(l_param)))
prop_P2=np.zeros((len(l_param)))
    
i=0

fig, axs = plt.subplots(4, 1, figsize=(8, 10))

# Psupply range loop
for param in l_param:
    if grazing == 'nograzing':
        [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply_arr*param,time, dt, gmax1=0,gmax2=0)
    if grazing == 'equalgrazing':
        [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply_arr*param,time, dt, gmax2=3.89,kZ2=5)
    if grazing == 'diffgrazing':
        [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply_arr*param,time, dt)
        
    # Record the equilibrium value of each state variable
    P_1[i] = P1[len(P1)-1]
    P_2[i] = P2[len(P2)-1]
    ratio[i] = P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1])
    P_O4[i] = PO4[len(PO4)-1]
    Z_[i] = Z[len(Z)-1]
    i=i+1
    print(i)
    
    # Check graphicly if the model is in steady state for each value of Psupply
    axs[0].plot(time, P1, label=f'P1 (param={param})', color="chartreuse")
    axs[1].plot(time, P2, label=f'P2 (param={param})',  color="green")
    axs[2].plot(time, Z, label=f'Z (param={param})', color="aqua")
    axs[3].plot(time, PO4, label=f'PO4 (param={param})', color="magenta")

plt.show()

filename_P1 = filename.replace('.txt', '_P1.txt')
filename_P2 = filename.replace('.txt', '_P2.txt')
filename_ratio = filename.replace('.txt', '_ratio.txt')
filename_PO4 = filename.replace('.txt', '_PO4.txt')
filename_Z = filename.replace('.txt', '_Z.txt')

np.savetxt(filename_P1, P_1)
np.savetxt(filename_P2, P_2)
np.savetxt(filename_ratio, ratio)
np.savetxt(filename_PO4, P_O4)
np.savetxt(filename_Z, Z_)