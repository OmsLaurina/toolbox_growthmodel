"""
Plot the steady-state outputs deviation from default values as a function of the value chosen for Psupply. Create the figure 7 of the paper.

"""

import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science','no-latex'])
plt.close('all')

# Load the std_variations_abs data from the outputs_deviation.py code
test_values = ['0.01', '0.1']
variations = []
for test in test_values:
    variations.append(np.loadtxt(f'sensitivity_variations_{test}.txt'))

plt.rc('font', size=10)
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
plt.subplots_adjust(wspace=0, hspace=0)

param_names = [r'$\mu^{max}_1$', r'$\mu^{max}_2$', r'$K^P_1$', r'$K^P_2$', r'$g^{max}_1$', r'$g^{max}_2$', r'$K^{Z}_1$', r'$K^{Z}_2$']

# Loop on test values
for i, test in enumerate(test_values):
    variations_abs = variations[i]
    cax = axs[i].matshow(variations_abs, cmap='Spectral')
    axs[i].set_xticks(np.arange(4))
    axs[i].set_yticks(np.arange(len(param_names)))
    axs[i].set_xticklabels([r'$PO4$', r'$P_1$', r'$P_2$', r'$R$'])
    axs[i].set_yticklabels(param_names)
    axs[i].xaxis.tick_bottom()
    axs[i].set_xlabel('State variable')
    axs[i].set_ylabel('Parameter')
    axs[i].set_title(rf'$\mathbf{{P_{{supply}} = {test}}}$', fontweight='bold', fontsize=12)
    
    if i == 0:
        cbar = plt.colorbar(cax, ax=axs[i], 
                            ticks=[-100, 0, 100], format='%.0f%%', extend='both')
        cax.set_clim(-100, 100)
    elif i == 1:
        cbar = plt.colorbar(cax, ax=axs[i], 
                            ticks=[-900, 0, 900], format='%.0f%%', extend='both')
        cax.set_clim(-900, 900)
    
    cbar.set_label('Absolute variation [%]', fontsize=10)
    
plt.tight_layout()
plt.savefig("../figures/deviation.pdf", format='pdf')
plt.show()