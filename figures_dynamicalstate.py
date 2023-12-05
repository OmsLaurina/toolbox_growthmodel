"""
* Plot the dynamical-state outputs of the model as a function of the value chosen for the amplitude b.
* The configuration is defined within the "set_up" function in the "set_up.py" script.
* First you need to create the outputs from the "outputs_dynamicalstate.py" script.
"""

import matplotlib.pyplot as plt
import numpy as np
import scienceplots
import sys
sys.path.append('../')

plt.style.use(['science','no-latex'])
plt.close('all')

# 1 pulse
# Load data created from outputs_dynamicalstate.py
data1 = np.loadtxt('dynamicalstate_pulse1_b0.08_d90.txt')
time_flux, P1, P2, Z, PO4, Psupply = data1.T

plt.figure(1)
name_pdf = '../figures/dynamicalstate_temporalevolution_pulse1_b0.08_d90.pdf'
plt.rc('font', size=7) 
plt.plot(time_flux,Psupply, label=r'$P_{supply}$', color="darkgray")
plt.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time_flux, P2, label=r'$P_2$', color="green")
plt.plot(time_flux, Z, label= r'$Z$', color="aqua")
plt.plot(time_flux,PO4, label=r'$PO4$', color="magenta")
plt.xlabel('Time [day]',fontsize=10)
plt.ylabel('Masses  [mmolC.m\u207B\u00B3]',fontsize=10)
plt.ylim(0,0.6)
plt.legend(frameon=True)
legend = plt.legend(frameon=True, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title(' "One pulses" ',fontweight='bold',fontsize=12)
plt.savefig(name_pdf, format='pdf')

# 2 pulses
# Load data created from outputs_dynamicalstate.py
data2 = np.loadtxt('dynamicalstate_pulse2_b0.08_d90.txt')
time_flux, P1, P2, Z, PO4, Psupply = data2.T

plt.figure(2)
name_pdf = '../figures/dynamicalstate_temporalevolution_pulse2_b0.08_d90.pdf'
plt.rc('font', size=7) 
plt.plot(time_flux,Psupply, label=r'$P_{supply}$', color="darkgray")
plt.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time_flux, P2, label=r'$P_2$', color="green")
plt.plot(time_flux, Z, label= r'$Z$', color="aqua")
plt.plot(time_flux,PO4, label=r'$PO4$', color="magenta")
plt.xlabel('Time [day]',fontsize=10)
plt.ylabel('Masses  [mmolC.m\u207B\u00B3]',fontsize=10)
plt.ylim(0,0.6)
plt.legend(frameon=True)
legend = plt.legend(frameon=True, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title(' "Two pulses" ',fontweight='bold',fontsize=12)
plt.savefig(name_pdf, format='pdf')

# 3 pulses
plt.figure(3)
# Load data created from outputs_dynamicalstate.py
data3 = np.loadtxt('dynamicalstate_pulse3_b0.08_d90.txt')
time_flux, P1, P2, Z, PO4, Psupply = data3.T

name_pdf = '../figures/dynamicalstate_temporalevolution_pulse3_b0.08_d90.pdf'
plt.rc('font', size=7) 
plt.plot(time_flux,Psupply, label=r'$P_{supply}$', color="darkgray")
plt.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time_flux, P2, label=r'$P_2$', color="green")
plt.plot(time_flux, Z, label= r'$Z$', color="aqua")
plt.plot(time_flux,PO4, label=r'$PO4$', color="magenta")
plt.xlabel('Time [day]',fontsize=10)
plt.ylabel('Masses  [mmolC.m\u207B\u00B3]',fontsize=10)
plt.ylim(0,0.6)
plt.legend(frameon=True)
legend = plt.legend(frameon=True, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title(' "Three pulse" ',fontweight='bold',fontsize=12)
plt.savefig(name_pdf, format='pdf')