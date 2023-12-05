"""
* Calculate and record the steady-state outputs of the model as a function of a range of Psupply value.
* The configuration is defined within the "set_up" function in the "set_up.py" script.
"""

import numpy as npy
import numpy as np
from growth_model import growth_model
import matplotlib.pyplot as plt
import scienceplots
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from set_up import set_up_steadystate

plt.style.use(['science', 'no-latex'])
plt.close('all')

# Configuration
dt, end_time, time, Psupply, Psupply_cst, Psupply_senstest, min_param, max_param,l_param = set_up_steadystate()

# Choose the grazing control type
grazing = 'diffgrazing'

# Parameter and file names
param = 'Psupp'
name_param = r'$P_{\mathrm{supply}}$ [mmolC m$^{-3}d^{-1}$]'
filename = f"../outputs/steadystate_senstitivitytest_{param}_{grazing}.txt"

# Create the empty matrices
P_1=np.zeros((len(l_param)))
P_2=np.zeros((len(l_param)))
ratio=np.zeros((len(l_param)))
P_O4=np.zeros((len(l_param)))
Z_=np.zeros((len(l_param)))

i=0

fig, axs = plt.subplots(2, 2, figsize=(8, 6))
fig2 = plt.figure(figsize=(8, 6))
ax2 = fig2.add_subplot(111, projection='3d')

# Psupply range loop
for param in l_param:
    if grazing == 'nograzing':
        [P1,P2,Z,PO4,arg]=growth_model(Psupply_senstest*param,time, dt, gmax1=0,gmax2=0)
    if grazing == 'equalgrazing':
        [P1,P2,Z,PO4,arg]=growth_model(Psupply_senstest*param,time, dt, gmax2=3.89,kZ2=5)
    if grazing == 'diffgrazing':
        [P1,P2,Z,PO4,arg]=growth_model(Psupply_senstest*param,time, dt)
        
    # Record the equilibrium value of each state variable
    P_1[i] = P1[len(P1)-1]
    P_2[i] = P2[len(P2)-1]
    ratio[i] = P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1])
    P_O4[i] = PO4[len(PO4)-1]
    Z_[i] = Z[len(Z)-1]
    
    i=i+1
    print(i)
    
    # Check graphicly if the model is in steady state for each value of Psupply
    
    color = plt.cm.plasma(i / len(l_param))
    titles = [r'$PO_4$', r'$P_1$',r'$P_2$',r'$Z$']
    # Temporal evolution
    axs[0,0].plot(time, PO4, label=f'PO4 (param={param})', color=color)
    axs[0,1].plot(time, P1, label=f'P1 (param={param})', color=color)
    axs[1,0].plot(time, P2, label=f'P2 (param={param})', color=color)
    axs[1,1].plot(time, Z, label=f'Z (param={param})', color=color)

    # Trajectories in 3D
    ax2.plot(P1, P2, Z, label=f'{name_param} = {param:.2f}', color=color)
    
sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma)
sm.set_array(l_param)

cbar = fig.colorbar(sm, ax=axs, orientation='vertical', shrink=0.5, pad=0.05)
cbar.set_label(name_param)

for i in range(2):
    for j in range(2):
        axs[i, j].set_title(titles[i * 2 + j])
for i in range(2):
    axs[1, i].set_xlabel("Time [d]")
for i in range(2):
    axs[i, 0].set_ylabel("Masses [mmolC m$^{-3}$]")


fig.savefig('../figures/evol_temp.pdf', format='pdf')
plt.show()

ax2.set_xlabel(r'$P_1$ [mmolC m$^{-3}$]')
ax2.set_ylabel(r'$P_2$ [mmolC m$^{-3}$]')
ax2.set_zlabel(r'$Z$ [mmolC m$^{-3}$]')
ax2.scatter(P1[0], P2[0], Z[0], marker = '*', c='k')
sm = ScalarMappable(cmap=plt.cm.plasma)
sm.set_array(l_param)
cbar = plt.colorbar(sm, orientation='horizontal', shrink=0.5, pad=0.05)
cbar.set_label(name_param)
fig2.savefig('../figures/trajectories.pdf', format='pdf')
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