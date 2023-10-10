"""

Plot the bifurcation diagrams in function of the psupply range.

"""
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
import sys
from growth_model_2P1Z_v10_simplify import growth_model_2P1Z_v10_simplify

sys.path.append('../')
plt.style.use(['science','no-latex'])
plt.close('all')


l_param = np.loadtxt('../outputs/l_param.txt')
name_param = r'$P_{SUPPLY}$ $[mmolCm^{-3}d^{-1}]$'
plt.figure(figsize=(7.5, 3))

# Left panel
# Load data from outputs_jacobianmatrix_simplifymodel.py code
equilibrium = ['P1', 'P2']
real_parts_list = []
max_real_parts_list = []

for eq in equilibrium:
    real_parts_list.append(np.loadtxt(f'../outputs/real_parts{eq}.txt'))
    max_real_parts_list.append(np.loadtxt(f'../outputs/max_real_parts{eq}.txt'))

plt.subplot(1, 2, 1)
for i, eq in enumerate(equilibrium):
    real_parts = real_parts_list[i]
    max_real_parts = max_real_parts_list[i]

    boolean_values = np.all(real_parts < 0, axis=1)
    stable_color = 'chartreuse' if i == 0 else 'darkgreen'
    instable_color = 'palegreen' if i == 0 else 'darkseagreen'
    plt.axhline(y=0, color='red', linestyle='--')
    plt.scatter(l_param, max_real_parts, c=np.where(boolean_values, stable_color, instable_color), s=5)
    plt.xlabel(name_param)
    plt.ylabel(r'$\lambda^{max}$')

# Right panel
# Load data from outputs_jacobianmatrix_fullmodel.py code
real_parts2 = np.loadtxt('../outputs/real_parts2.txt')
max_real_parts2 = np.loadtxt('../outputs/max_real_parts2.txt')

plt.subplot(1, 2, 2)
boolean_values = np.all(real_parts2 < 0, axis=1)
stable_color = 'black'
instable_color = 'grey'
plt.axhline(y=0, color='red', linestyle='--')
plt.plot(l_param, max_real_parts2, 'lightgray')
plt.scatter(l_param, max_real_parts2, c=np.where(boolean_values, stable_color, instable_color), s=5)
plt.xlabel(name_param)
plt.ylabel(r'$\lambda^{max}$')

plt.tight_layout()
plt.savefig('../figures/diagrbifurc.pdf', format='pdf')
plt.show()