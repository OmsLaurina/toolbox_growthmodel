#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 09:36:07 2023

@author: loms
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import math
import scienceplots
import sys
from growth_model import growth_model
from growth_model_RK45 import growth_model

sys.path.append('../')
plt.style.use(['science','no-latex'])
plt.close('all')

# Define the set up
dt = 0.1
end_time = 2000
time = np.arange(0,end_time,dt)
Psupply_moy = 1
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

# Number of value tested
n = 100  
param = 'psupp'
# name for plot
name_param= r'$P_{\mathrm{supply}}$ [mmolCm$^{-3}d^{-1}$]'

#Tested range of values
min_param = 0.01
max_param = 0.4
l_param = np.linspace(min_param,max_param,n)

eigenvalues=[]
real_parts_list=[]
max_real_parts = []

i=0

equilibre = 'P2null'
for param in l_param:
                                    
    a = arg['gamma']*(1-arg['epsilon2Z']*(1-arg['gamma']))
    b = -(1-arg['epsilon2Z']*(1-arg['gamma'])+arg['epsilon1Z']*arg['gamma'])*arg['m1Z']
    c = arg['epsilon1Z']*arg['m1Z']**2-Psupply_moy*param*arg['m2Z']
    
    delta = b**2-4*a*c
    x1 = (-b+math.sqrt(delta))/(2*a)
    x2 =(-b-math.sqrt(delta))/(2*a)
    
    if equilibre == 'P2null':
        
        name_pdf_real_parts = '../outputs/real_partsP1_full.txt'
        name_pdf_max_real_parts = '../outputs/max_real_partsP1_full.txt'
        
        P1_barre = x1*arg['kZ1']/(arg['gmax1']-x1)
        P2_barre = 0
        Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
        PO4_barre = arg['kP1']*(x1*Z_barre+arg['mP1']*P1_barre)/(arg['umax1']*P1_barre-(x1*Z_barre+arg['mP1']*P1_barre))
        
        if P1_barre > 0 and Z_barre > 0 and PO4_barre > 0:
            print(f"Équilibre 'P2null' est positif pour param={param}")
        else:
            print(f"Équilibre 'P2null' n'est pas positif pour param={param}")
    
    if equilibre == 'P1null':
        
        name_pdf_real_parts = '../outputs/real_partsP2_full.txt'
        name_pdf_max_real_parts = '../outputs/max_real_partsP2_full.txt'
        
        P1_barre = 0
        P2_barre = x1*arg['kZ2']/(arg['gmax2']-x1)
        Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
        PO4_barre = arg['kP2']*(x1*Z_barre+arg['mP2']*P2_barre)/(arg['umax2']*P2_barre-(x1*Z_barre+arg['mP2']*P2_barre))
        
        if P2_barre > 0 and Z_barre > 0 and PO4_barre > 0:
            print(f"Équilibre {equilibre} est positif pour param={param}")
        else:
            print(f"Équilibre {equilibre} n'est pas positif pour param={param}")
        
    if equilibre == 'coexistence':
        
        name_pdf_real_parts = '../outputs/real_partscoex_full.txt'
        name_pdf_max_real_parts = '../outputs/max_real_partscoex_full.txt'
        
        [P1,P2,Z,PO4,arg]=growth_model(Psupply_arr*param, time, dt)
        P1_barre = P1[-1]
        P2_barre = P2[-1]
        if P2_barre<0.02:
            P2_barre = 0
        Z_barre = Z[-1]
        PO4_barre = PO4[-1]
        
        if P1_barre > 0 and P2_barre > 0 and Z_barre > 0 and PO4_barre > 0:
            print(f"Équilibre {equilibre} est positif pour param={param}")
        else:
            print(f"Équilibre {equilibre} n'est pas positif pour param={param}")
            
    
    # Jacobian matrix
    
    j11 = - P1_barre * arg['umax1'] * arg['kP1'] / (PO4_barre + arg['kP1'])**2 - P2_barre * arg['umax2'] * arg['kP2'] / (PO4_barre + arg['kP2'])**2
    j12 = Z_barre *arg['epsilon2Z']* (1 -arg['gamma']) * (arg['gmax1'] * (P2_barre + arg['kZ1']) / (P1_barre + P2_barre + arg['kZ1'])**2 - arg['gmax2']* P2_barre / (P1_barre + P2_barre + arg['kZ2'])**2) + arg['mP1'] - arg['umax1'] * PO4_barre / (PO4_barre + arg['kP1'])
    j13 = Z_barre *arg['epsilon2Z']* (1 -arg['gamma']) * (arg['gmax2']* (P1_barre + arg['kZ2']) / (P1_barre + P2_barre + arg['kZ2'])**2 - arg['gmax1'] * P1_barre / (P1_barre + P2_barre + arg['kZ1'])**2) + arg['mP2']- arg['umax2'] * PO4_barre / (PO4_barre + arg['kP2'])
    j14 = arg['epsilon2Z']* (1 -arg['gamma']) * (arg['gmax1'] * P1_barre / (P1_barre + P2_barre + arg['kZ1']) + arg['gmax2']* P2_barre / (P1_barre + P2_barre + arg['kZ2'])) + arg['epsilon1Z']*arg['m1Z']
    
    j21 = (P1_barre * arg['umax1'] * arg['kP1']) / (PO4_barre + arg['kP1'])**2
    j22 = (arg['umax1'] * PO4_barre) / (PO4_barre + arg['kP1']) - (Z_barre * arg['gmax1'] * (P2_barre + arg['kZ1'])) / ((P1_barre + P2_barre + arg['kZ1'])**2) - arg['mP1']
    j23 = (Z_barre * arg['gmax1'] * P1_barre) / (P1_barre + P2_barre + arg['kZ1'])**2
    j24 = - (arg['gmax1'] * P1_barre) / (P1_barre + P2_barre + arg['kZ1'])

    j31 = (P2_barre * arg['umax2'] * arg['kP2']) / (PO4_barre + arg['kP2'])**2
    j32 = (Z_barre *arg['gmax2']* P2_barre) / (P1_barre + P2_barre + arg['kZ1'])**2
    j33 = (arg['umax2'] * PO4_barre) / (PO4_barre + arg['kP2']) - (Z_barre *arg['gmax2']* (P1_barre + arg['kZ2'])) / (P1_barre + P2_barre + arg['kZ1'])**2 - arg['mP2']
    j34 = - (arg['gmax2']* P2_barre) / (P1_barre + P2_barre + arg['kZ2'])

    j41 = 0
    j42 = arg['gamma'] * Z_barre * ((arg['gmax1'] * (P2_barre + arg['kZ1'])) / (P1_barre + P2_barre + arg['kZ1'])**2 - (arg['gmax2']* P2_barre) / (P1_barre + P2_barre + arg['kZ2'])**2)
    j43 = arg['gamma'] * Z_barre * ((arg['gmax2']* (P1_barre + arg['kZ2'])) / (P1_barre + P2_barre + arg['kZ2'])**2 - (arg['gmax1'] * P1_barre) / (P1_barre + P2_barre + arg['kZ1'])**2)
    j44 = arg['gamma'] * ((arg['gmax1'] * P1_barre) / (P1_barre + P2_barre + arg['kZ1']) + (arg['gmax2']* P2_barre) / (P1_barre + P2_barre + arg['kZ2'])) - 2*arg['m2Z']*Z_barre-arg['m1Z']
    
    J = [[j11,j12,j13,j14], [j21,j22,j23,j24], [j31,j32,j33,j34],[j41,j42,j43,j44]]
    
    eigenvalues.append(np.linalg.eigvals(J))  

    real_parts = np.real(eigenvalues[-1])
    real_parts_list.append(real_parts) 
    max_real_part = np.max(real_parts)
    max_real_parts.append(max_real_part)
    
    if np.all(real_parts < 0):
        print("Le modèle est stable")
    else:
        print("Le modèle est instable")
    real_parts_matrix = np.vstack(real_parts_list)
    
    i=i+1
    print(i)
 
np.savetxt(name_pdf_real_parts, real_parts_matrix)
np.savetxt(name_pdf_max_real_parts, max_real_parts)

# Check if the solutions obtain with the jacobian matrix are in consistence with numerical solutions

if equilibre == 'P2null':
    
    p=0.01
    Psupply = [p] * len(time)
    [P1,P2,Z,PO4,arg]=growth_model(Psupply, time, dt, P1_ini=0.6, P2_ini=0.1)
    
    a = arg['gamma']*(1-arg['epsilon2Z']*(1-arg['gamma']))
    b = -(1-arg['epsilon2Z']*(1-arg['gamma'])+arg['epsilon1Z']*arg['gamma'])*arg['m1Z']
    c = arg['epsilon1Z']*arg['m1Z']**2-p*arg['m2Z']

    delta = b**2-4*a*c
    x1 = (-b+math.sqrt(delta))/(2*a)
    x2 =(-b-math.sqrt(delta))/(2*a)
    
    P1_barre = x1*arg['kZ1']/(arg['gmax1']-x1)
    P2_barre = 0
    Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
    PO4_barre = arg['kP1']*(x1*Z_barre+arg['mP1']*P1_barre)/(arg['umax1']*P1_barre-(x1*Z_barre+arg['mP1']*P1_barre))

if equilibre == 'P1null':
    
    p=0.3
    Psupply = [p] * len(time)
    [P1,P2,Z,PO4,arg]=growth_model(Psupply, time, dt, P1_ini=0.6, P2_ini=0.1)
    
    a = arg['gamma']*(1-arg['epsilon2Z']*(1-arg['gamma']))
    b = -(1-arg['epsilon2Z']*(1-arg['gamma'])+arg['epsilon1Z']*arg['gamma'])*arg['m1Z']
    c = arg['epsilon1Z']*arg['m1Z']**2-p*arg['m2Z']

    delta = b**2-4*a*c
    x1 = (-b+math.sqrt(delta))/(2*a)
    
    P1_barre = 0
    P2_barre = x1*arg['kZ2']/(arg['gmax2']-x1)
    Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
    PO4_barre = arg['kP2']*(x1*Z_barre+arg['mP2']*P2_barre)/(arg['umax2']*P2_barre-(x1*Z_barre+arg['mP2']*P2_barre))

if equilibre == 'coexistence':
    
    p=0.05
    Psupply = [p] * len(time)
    
    [P1,P2,Z,PO4,arg]=growth_model(Psupply, time, dt)
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


# Check if the solutions obtain with the jacobian matrix are in consistence with numerical solutions from RK45

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
if equilibre == 'P2null':
    P1_barre = x1*arg['kZ1']/(arg['gmax1']-x1)
    P2_barre = 0
    Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
    PO4_barre = arg['kP1']*(x1*Z_barre+arg['mP1']*P1_barre)/(arg['umax1']*P1_barre-(x1*Z_barre+arg['mP1']*P1_barre))

if equilibre == 'P1null':
    P1_barre = 0
    P2_barre = x1*arg['kZ2']/(arg['gmax2']-x1)
    Z_barre = (arg['gamma']*x1-arg['m1Z'])/arg['m2Z']
    PO4_barre = arg['kP2']*(x1*Z_barre+arg['mP2']*P2_barre)/(arg['umax2']*P2_barre-(x1*Z_barre+arg['mP2']*P2_barre))
    
if equilibre == 'coexistence':
    [P1,P2,Z,PO4,arg]=growth_model(Psupply, time,  dt)
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
