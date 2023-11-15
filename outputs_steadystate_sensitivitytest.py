"""

Calculate and record the steady-state outputs of the model as a function of a range of Psupply value.

"""
import numpy as npy
import numpy as np
from growth_model import growth_model
import matplotlib.pyplot as plt
import scienceplots
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
name_param = r'$P_{\mathrm{supply}}$ [mmolC m$^{-3}d^{-1}$]'

# Tested range of values
min_param = 0.01
max_param = 0.1
l_param = np.linspace(min_param,max_param,n)
    
filename = f"../outputs/steadystate_senstitivitytest_{param}_{grazing}.txt"

# Create the empty matrices
P_1=np.zeros((len(l_param)))
P_2=np.zeros((len(l_param)))
ratio=np.zeros((len(l_param)))
P_O4=np.zeros((len(l_param)))
Z_=np.zeros((len(l_param)))

i=0

fig, axs = plt.subplots(2, 2, figsize=(6, 4))

fig2 = plt.figure(figsize=(8, 6))
ax2 = fig2.add_subplot(111, projection='3d')

# Psupply range loop
for param in l_param:
    if grazing == 'nograzing':
        [P1,P2,Z,PO4,arg]=growth_model(Psupply_arr*param,time, dt, gmax1=0,gmax2=0)
    if grazing == 'equalgrazing':
        [P1,P2,Z,PO4,arg]=growth_model(Psupply_arr*param,time, dt, gmax2=3.89,kZ2=5)
    if grazing == 'diffgrazing':
        [P1,P2,Z,PO4,arg]=growth_model(Psupply_arr*param,time, dt)
        
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
    titles = [r'$PO4$', r'$P_1$',r'$P_2$',r'$Z$']
    # Temporal evolution
    axs[0,0].plot(time, PO4, label=f'PO4 (param={param})', color=color)
    axs[0,1].plot(time, P1, label=f'P1 (param={param})', color=color)
    axs[1,0].plot(time, P2, label=f'P2 (param={param})', color=color)
    axs[1,1].plot(time, Z, label=f'Z (param={param})', color=color)

    # Trajectories in 3D
    ax2.plot(P1, P2, Z, label=f'{name_param} = {param:.2f}', color=color)
    
sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma)
sm.set_array(l_param)

# Créez une seule colorbar pour l'ensemble du subplot
cbar = fig.colorbar(sm, ax=axs, orientation='vertical', shrink=0.5, pad=-0.45)
cbar.set_label(name_param)

for i in range(2):
    for j in range(2):
        axs[i, j].set_title(titles[i * 2 + j])
plt.tight_layout()

for i in range(2):
    axs[1, i].set_xlabel("Time [d]")
for i in range(2):
    axs[i, 0].set_ylabel("Masses [mmolC m$^{-3}$]")
plt.subplots_adjust(wspace=0.3)
plt.show()
plt.savefig('../figures/evol_temp.pdf', format='pdf')

ax2.set_xlabel(r'$P_1$ [mmolC m$^{-3}$]')
ax2.set_ylabel(r'$P_2$ [mmolC m$^{-3}$]')
ax2.set_zlabel(r'$Z$ [mmolC m$^{-3}$]')
ax2.scatter(P1[0], P2[0], Z[0], marker = '*', c='k')
sm = ScalarMappable(cmap=plt.cm.plasma)
sm.set_array(l_param)
cbar = plt.colorbar(sm, orientation='horizontal', shrink=0.5, pad=0.05)
cbar.set_label(name_param)
plt.show()
plt.savefig('../figures/trajectories.pdf', format='pdf')

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