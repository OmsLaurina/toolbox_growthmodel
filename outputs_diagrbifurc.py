#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 09:36:07 2023

@author: loms
"""

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as npy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import math
import scienceplots
import sys
from growth_model_2P1Z_v10 import growth_model_2P1Z_v10 
from growth_model_2P1Z_v10_RK45 import growth_model
sys.path.append('../')
plt.style.use(['science','no-latex'])
plt.close('all')


# Timestep (same than timestep of the growth function : by default is 0.2)
dt = 0.1
# Choose the time at which the simulation end (days)
end_time = 2000
# create a time vector from 0 to end_time, seperated by a timestep
time = npy.arange(0,end_time,dt)
Psupply_moy = 1
p=0.01
Psupply = [p] * len(time)
Psupply_arr = np.array([Psupply_moy]*len(time))


# Parameters
arg = {
    'umax1':1.9872,   # maximum growth rates of P1 [d^{-1}]
    'umax2':2.7648,   # maximum growth rates of P2 [d^{-1}]
    'gmax1':3.89,     # maximum grazing rates of Z on P1 [d^{-1}]
    'gmax2':0.43,     # maximum grazing rates of Z on P2 [d^{-1}]
    'kP1':1,          # half-saturation constant for P1 on PO4 [mmolC m^{-3}]
    'kP2':3,	      # half-saturation constant for P2 on PO4 [mmolC m^{-3}]
    'kZ1':5,		  # half-saturation constant for Z on P1 [mmolC m^{-3}]
    'kZ2':20,		  # half-saturation constant for Z on P2 [mmolC m^{-3}]
    'mP1':0.1,	      # P1 natural mortality rate [d^{-1}]
    'mP2':0.2,	      # P2 natural mortality rate [d^{-1}]
    'm1Z':0.1,	      # Z natural mortality rate [d^{-1}]
    'm2Z':0.061,	  # Z quadratic mortality rate [mmolC^{-1} m^{3} d^{-1}]
    'gamma':0.6,      # conversion factor from P to Z [/]
    'epsilonP':1,     # fraction of P natural mortality that is available as regenerated PO4 [/]
    'epsilon1Z':0.3,  # fraction of Z natural mortality that is available as regenerated PO4 [/]
    'epsilon2Z':0.7,  # fraction of Z excretion that is available as regenerated PO4 [/]
    'P1_ini':0.6,     # initial biomass of P1 [mmolC m^{-3}]
    'P2_ini':0.1,     # initial biomass of P2 [mmolC m^{-3}]
    'Z_ini':0.6,      # initial biomass of Z [mmolC m^{-3}]
    'PO4_ini':0.5,    # initial biomass of PO4 [mmolC m^{-3}]
    }

[P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply, time, dt,P1_ini=0.6, P2_ini=0.1)

print(P1[-1])
print(P2[-1])

n = 100 #number of value tested
#name of the parameters in the model you want to study 
param1 = 'psupp'
# name for plot
name_param1= r'$P_{\mathrm{supply}}$ [mmolCm$^{-3}d^{-1}$]'

#Tested range of values
min_param1 = 0.01
max_param1 = 0.12
l_param1 = np.linspace(min_param1,max_param1,n)

eigenvalues=[]
real_parts_list=[]
max_real_parts = []

i=0

for param1 in l_param1:
                                    
    a = arg['gamma']*(1-arg['epsilon2Z']*(1-arg['gamma']))
    b = -(1-arg['epsilon2Z']*(1-arg['gamma'])+arg['epsilon1Z']*arg['gamma'])*arg['m1Z']
    c = arg['epsilon1Z']*arg['m1Z']**2-Psupply_moy*param1*arg['m2Z']
    
    delta = b**2-4*a*c
    x1 = (-b+math.sqrt(delta))/(2*a)
    x2 =(-b-math.sqrt(delta))/(2*a)
    
    if arg['P2_ini'] == 0: 
        P1_barre = x1*arg['kZ1']/(arg['gmax1']-x1)
        P2_barre = 0
        Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
        PO4_barre = arg['kP1']*(x1*Z_barre+arg['mP1']*P1_barre)/(arg['umax1']*P1_barre-(x1*Z_barre+arg['mP1']*P1_barre))
    
    if arg['P1_ini'] == 0: 
        P1_barre = 0
        P2_barre = x1*arg['kZ2']/(arg['gmax2']-x1)
        Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
        PO4_barre = arg['kP2']*(x1*Z_barre+arg['mP2']*P2_barre)/(arg['umax2']*P2_barre-(x1*Z_barre+arg['mP2']*P2_barre))
        
    if arg['P1_ini'] == 0.6 and arg['P2_ini'] == 0.1:
        [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply_arr*param1, time, dt,P1_ini=0.6, P2_ini=0.1)
        P1_barre = P1[-1]
        P2_barre = P2[-1]
        Z_barre = Z[-1]
        PO4_barre = PO4[-1]
    
    # Matrice Jacobienne
    
    j11 = (arg['umax1'] * PO4_barre) / (PO4_barre + arg['kP1']) - (Z_barre * arg['gmax1'] * (P2_barre + arg['kZ1'])) / ((P1_barre + P2_barre + arg['kZ1'])**2) - arg['mP1']
    j12 = (Z_barre * arg['gmax1'] * P1_barre) / (P1_barre + P2_barre + arg['kZ1'])**2
    j13 = - (arg['gmax1'] * P1_barre) / (P1_barre + P2_barre + arg['kZ1'])
    j14 = (P1_barre * arg['umax1'] * arg['kP1']) / (PO4_barre + arg['kP1'])**2

    j21 = (Z_barre *arg['gmax2']* P2_barre) / (P1_barre + P2_barre + arg['kZ1'])**2
    j22 = (arg['umax2'] * PO4_barre) / (PO4_barre + arg['kP2']) - (Z_barre *arg['gmax2']* (P1_barre + arg['kZ1'])) / (P1_barre + P2_barre + arg['kZ1'])**2 - arg['mP2']
    j23 = - (arg['gmax2']* P2_barre) / (P1_barre + P2_barre + arg['kZ1'])
    j24 = (P2_barre * arg['umax2'] * arg['kP2']) / (PO4_barre + arg['kP2'])**2

    j31 = arg['gamma'] * Z_barre * ((arg['gmax1'] * (P2_barre + arg['kZ1'])) / (P1_barre + P2_barre + arg['kZ1'])**2 - (arg['gmax2']* P2_barre) / (P1_barre + P2_barre + arg['kZ1'])**2)
    j32 = arg['gamma'] * Z_barre * ((arg['gmax2']* (P1_barre + arg['kZ1'])) / (P1_barre + P2_barre + arg['kZ1'])**2 - (arg['gmax1'] * P1_barre) / (P1_barre + P2_barre + arg['kZ1'])**2)
    j33 = arg['gamma'] * ((arg['gmax1'] * P1_barre) / (P1_barre + P2_barre + arg['kZ1']) + (arg['gmax2']* P2_barre) / (P1_barre + P2_barre + arg['kZ1'])) - 2*arg['m2Z']*Z_barre-arg['m1Z']
    j34 = 0

    j41 = Z_barre *arg['epsilon2Z']* (1 -arg['gamma']) * (arg['gmax1'] * (P2_barre + arg['kZ1']) / (P1_barre + P2_barre + arg['kZ1'])**2 - arg['gmax2']* P2_barre / (P1_barre + P2_barre + arg['kZ1'])**2) + arg['mP1'] - arg['umax1'] * PO4_barre / (PO4_barre + arg['kP1'])
    j42 = Z_barre *arg['epsilon2Z']* (1 -arg['gamma']) * (arg['gmax2']* (P1_barre + arg['kZ1']) / (P1_barre + P2_barre + arg['kZ1'])**2 - arg['gmax1'] * P1_barre / (P1_barre + P2_barre + arg['kZ1'])**2) + arg['mP2']- arg['umax2'] * PO4_barre / (PO4_barre + arg['kP2'])
    j43 = arg['epsilon2Z']* (1 -arg['gamma']) * (arg['gmax1'] * P1_barre / (P1_barre + P2_barre + arg['kZ1']) + arg['gmax2']* P2_barre / (P1_barre + P2_barre + arg['kZ1'])) + arg['epsilon1Z']*arg['m1Z']
    j44 = - P1_barre * arg['umax1'] * arg['kP1'] / (PO4_barre + arg['kP1'])**2 - P2_barre * arg['umax2'] * arg['kP2'] / (PO4_barre + arg['kP2'])**2
    
    J = [[j11,j12,j13,j14], [j21,j22,j23,j24], [j31,j32,j33,j34],[j41,j42,j43,j44]]
    
    eigenvalues.append(np.linalg.eigvals(J))  

    real_parts = np.real(eigenvalues[-1])
    real_parts_list.append(real_parts) 
    max_real_part = np.max(real_parts)
    max_real_parts.append(max_real_part)
    t_res = 1/np.abs(max_real_part)
    
    if np.all(real_parts < 0): #un modèle est stable si toutes les valeurs propres de la MJ sont de parties réelles < 0 
        print("Le modèle est stable")
    else:
        print("Le modèle est instable")
    # print("Valeurs propres :", eigenvalues)
    # print("Parties réelles :", real_parts)
    real_parts_matrix = np.vstack(real_parts_list)
    # print("Matrice des parties réelles des valeurs propres :", real_parts_matrix)
    
    i=i+1
    print(i)
    print(param1)
    # print(max_real_part)
    print(t_res)

""" Figures """

plt.figure(2)
boolean_values = np.all(real_parts_matrix < 0, axis=1)
print(boolean_values)
stable_color = 'black'
instable_color = 'grey'
plt.axhline(y=0, color='red', linestyle='--')
plt.plot(l_param1, max_real_parts,'lightgray')
plt.scatter(l_param1, max_real_parts, c=np.where(boolean_values, stable_color, instable_color), s=5)
plt.xlabel(name_param1)
plt.ylabel(r'$\lambda^{max}$')
plt.savefig('diagrbifurc_complexemodel.pdf', format='pdf')

#### Test si les solutions analytiques collent avec les solutions numériques 

a = arg['gamma']*(1-arg['epsilon2Z']*(1-arg['gamma']))
b = -(1-arg['epsilon2Z']*(1-arg['gamma'])+arg['epsilon1Z']*arg['gamma'])*arg['m1Z']
c = arg['epsilon1Z']*arg['m1Z']**2-p*arg['m2Z']

delta = b**2-4*a*c
x1 = (-b+math.sqrt(delta))/(2*a)
x2 =(-b-math.sqrt(delta))/(2*a)

# Calcul des points d'équilibre
if arg['P2_ini'] == 0: 
    P1_barre = x1*arg['kZ1']/(arg['gmax1']-x1)
    P2_barre = 0
    Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
    PO4_barre = arg['kP1']*(x1*Z_barre+arg['mP1']*P1_barre)/(arg['umax1']*P1_barre-(x1*Z_barre+arg['mP1']*P1_barre))

if arg['P1_ini'] == 0: 
    P1_barre = 0
    P2_barre = x1*arg['kZ2']/(arg['gmax2']-x1)
    Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
    PO4_barre = arg['kP2']*(x1*Z_barre+arg['mP2']*P2_barre)/(arg['umax2']*P2_barre-(x1*Z_barre+arg['mP2']*P2_barre))

if arg['P1_ini'] == 0.6 and arg['P2_ini'] == 0.1:
    [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply, time, dt, P1_ini=0.6, P2_ini=0.1)
    P1_barre = P1[-1]
    P2_barre = P2[-1]
    Z_barre = Z[-1]
    PO4_barre = PO4[-1]

plt.figure(1)
plt.rc('font', size=7)
plt.plot(time, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time, P2, label=r'$P_2$', color="green")
plt.plot(time, Z, label=r'$Z$', color="aqua")
plt.plot(time, PO4, label=r'$PO4$', color="magenta")
plt.axhline(y=P1_barre, color='chartreuse', linestyle='--')
plt.axhline(y=P2_barre, color='green', linestyle='--')
plt.axhline(y=PO4_barre, color='magenta', linestyle='--')
plt.axhline(y=Z_barre, color='aqua', linestyle='--')


#### Test si les solutions analytiques collent avec les solutions numériques issus de la methode RK45

P1_ini = arg['P1_ini']
P2_ini = arg['P2_ini']
Z_ini = arg['Z_ini']
PO4_ini = arg['PO4_ini']
CI = [P1_ini, P2_ini, Z_ini, PO4_ini]
t_span = (0, end_time)
solution = solve_ivp(growth_model, [0, end_time], CI, method='RK45', t_eval=time, args=(p,))
P1, P2, Z, PO4 = solution.y
a = arg['gamma']*(1-arg['epsilon2Z']*(1-arg['gamma']))
b = -(1-arg['epsilon2Z']*(1-arg['gamma'])+arg['epsilon1Z']*arg['gamma'])*arg['m1Z']
c = arg['epsilon1Z']*arg['m1Z']**2-p*arg['m2Z']
delta = b**2-4*a*c
x1 = (-b+math.sqrt(delta))/(2*a)
x2 =(-b-math.sqrt(delta))/(2*a)

#Calcul des points d'équilibre
if P2_ini == 0: 
    P1_barre = x1*arg['kZ1']/(arg['gmax1']-x1)
    P2_barre = 0
    Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
    PO4_barre = arg['kP1']*(x1*Z_barre+arg['mP1']*P1_barre)/(arg['umax1']*P1_barre-(x1*Z_barre+arg['mP1']*P1_barre))

if P1_ini == 0: 
    P1_barre = 0
    P2_barre = x1*arg['kZ2']/(arg['gmax2']-x1)
    Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
    PO4_barre = arg['kP2']*(x1*Z_barre+arg['mP2']*P2_barre)/(arg['umax2']*P2_barre-(x1*Z_barre+arg['mP2']*P2_barre))
    
if arg['P1_ini'] == 0.6 and arg['P2_ini'] == 0.1:
    [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply, time,  dt, P1_ini=0.6, P2_ini=0.1)
    P1_barre = P1[-1]
    P2_barre = P2[-1]
    Z_barre = Z[-1]
    PO4_barre = PO4[-1]

plt.figure(3)
plt.rc('font', size=7)
plt.plot(time, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time, P2, label=r'$P_2$', color="green")
plt.plot(time, Z, label=r'$Z$', color="aqua")
plt.plot(time, PO4, label=r'$PO4$', color="magenta")
plt.axhline(y=P1_barre, color='chartreuse', linestyle='--')
plt.axhline(y=P2_barre, color='green', linestyle='--')
plt.axhline(y=PO4_barre, color='magenta', linestyle='--')
plt.axhline(y=Z_barre, color='aqua', linestyle='--')
