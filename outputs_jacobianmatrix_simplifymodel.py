"""

Calculate and record the jacobien matrix terms of the simplify model, in function of the psupply range.

"""

import matplotlib.pyplot as plt
import numpy as npy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import math
import scienceplots
import sys
from growth_model_2P1Z_v10_simplify import growth_model_2P1Z_v10_simplify
sys.path.append('../')
plt.close('all')


Psupply_moy = 1

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

umax1 = arg['umax1']
umax2 = arg['umax2']
gmax1 = arg['gmax1']
gmax2 = arg['gmax2']
kP1 = arg['kP1']
kP2 = arg['kP2']
kZ1 = arg['kZ1']
kZ2 = arg['kZ2']
mP1 = arg['mP1']
mP2 = arg['mP2']
m1Z = arg['m1Z']
m2Z = arg['m2Z']
gamma = arg['gamma']
epsilonP = arg['epsilonP']
epsilon1Z = arg['epsilon1Z']
epsilon2Z = arg['epsilon2Z']

# Number of tested value 
n = 100
param = 'psupp'
 
#Tested range of values
min_param = 0.01
max_param = 0.12
l_param = np.linspace(min_param,max_param,n)

# Define the empty matrices    
eigenvalues=[]
real_parts_list=[]
max_real_parts = []
c_lim = []
    
# Set the equilibrium P1null or P2null
Phyto = 'P1null'

i=0

# Loop on param
for param in l_param:

    if Phyto == 'P2null':
        
        name_pdf_real_parts = '../outputs/real_partsP1.txt'
        name_pdf_max_real_parts = '../outputs/max_real_partsP1.txt'
        
        az = arg['gamma']*arg['gmax1']/(arg['m2Z']*arg['kZ1'])
        cz = arg['m1Z']/arg['m2Z']
        apo4 = arg['kP1']*arg['gamma']*arg['gmax1']**2/(arg['umax1']*arg['m2Z']*arg['kZ1']**2)
        cpo4 = arg['kP1']/arg['umax1']*(arg['mP1']-arg['m1Z']*arg['gmax1']/(arg['m2Z']*arg['kZ1']))
        a = arg['epsilon2Z']*(1-arg['gamma'])*arg['gmax1']/arg['kZ1']*az-apo4*arg['umax1']/arg['kP1']
        b = -arg['epsilon2Z']*(1-arg['gamma'])*arg['gmax1']/arg['kZ1']*cz+arg['mP1']-cpo4*arg['umax1']/arg['kP1']+arg['epsilon1Z']*arg['m1Z']*az
        c = Psupply_moy*param-arg['epsilon1Z']*arg['m1Z']*cz
        delta = b**2-4*a*c 
        x1 = (-b-math.sqrt(delta))/(2*a)
        x2 =(-b+math.sqrt(delta))/(2*a)
        
        #Equilibria points
        P1 = x1
        P2 = 0
        Z = az*P1-cz
        PO4 = apo4*P1+cpo4
        
    else:
        
        name_pdf_real_parts = '../outputs/real_partsP2.txt'
        name_pdf_max_real_parts = '../outputs/max_real_partsP2.txt'
        
        az = arg['gamma']*arg['gmax2']/(arg['m2Z']*arg['kZ2'])
        cz = arg['m1Z']/arg['m2Z']
        apo4 = arg['kP2']*arg['gamma']*arg['gmax2']**2/(arg['umax2']*arg['m2Z']*arg['kZ2']**2)
        cpo4 = arg['kP2']/arg['umax2']*(arg['mP2']-arg['m1Z']*arg['gmax2']/(arg['m2Z']*arg['kZ2']))
        a = arg['epsilon2Z']*(1-arg['gamma'])*arg['gmax2']/arg['kZ2']*az-apo4*arg['umax2']/arg['kP2']
        b = -arg['epsilon2Z']*(1-arg['gamma'])*arg['gmax2']/arg['kZ2']*cz+arg['mP2']-cpo4*arg['umax2']/arg['kP2']+arg['epsilon1Z']*arg['m1Z']*az
        c = Psupply_moy*param-arg['epsilon1Z']*arg['m1Z']*cz
        delta = b**2-4*a*c 
        x1 = (-b-math.sqrt(delta))/(2*a)
        x2 =(-b+math.sqrt(delta))/(2*a)
        
        #Equilibria points
        P1 = 0
        P2 = x1
        Z = az*P2-cz
        PO4 = apo4*P2+cpo4
    
    # Jacobian matrix
    j11 = (umax1 / kP1) * PO4 - (gmax1 / kZ1) * Z - mP1
    j12 = 0
    j13 = -(gmax1 / kZ1) * P1
    j14 = (umax1 / kP1) * P1
    j21 = 0
    j22 = (umax2 / kP2) * PO4 - (gmax2 / kZ2) * Z - mP2
    j23 = -(gmax2 / kZ2) * P2
    j24 = (umax2 / kP2) * P2
    j31 = Z * gamma * (gmax1 / kZ1)
    j32 = Z * gamma * (gmax2 / kZ2)
    j33 = gamma * ((gmax1 / kZ1) * P1 + (gmax2 / kZ2) * P2) - 2 * m2Z * Z - m1Z
    j34 = 0
    j41 = Z * epsilonP * (1 - gamma) * (gmax1 / kZ1) + mP1 - (umax1 / kP1) * PO4
    j42 = Z * epsilonP * (1 - gamma) * (gmax2 / kZ2) + mP2 - (umax2 / kP2) * PO4
    j43 = epsilonP * (1 - gamma) * ((gmax1 / kZ1) * P1 + (gmax2 / kZ2) * P2) + epsilon2Z * m2Z
    j44 = -(umax1 / kP1) * P1 - (umax2 / kP2) * P2
    J = [[j11,j12,j13,j14], [j21,j22,j23,j24], [j31,j32,j33,j34], [j41,j42,j43,j44]]
    
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
    # print("Matrice des parties réelles des valeurs propres :", real_parts_matrix)
    i=i+1

np.savetxt(name_pdf_real_parts, real_parts_matrix)
np.savetxt(name_pdf_max_real_parts, max_real_parts)
np.savetxt('../outputs/l_param.txt', l_param)

# Check if the solutions obtain with the jacobian matrix are in consistence with numerical solutions
# Define the set up
dt = 0.2
end_time = 2000
time = npy.arange(0,end_time,dt)


if Phyto == 'P2null':
    
    p = 0.01
    Psupply = [p] * len(time)
    
    [P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10_simplify(Psupply, time, P2_ini = 0)
       
    az = arg['gamma']*arg['gmax1']/(arg['m2Z']*arg['kZ1'])
    cz = arg['m1Z']/arg['m2Z']
    
    apo4 = arg['kP1']*arg['gamma']*arg['gmax1']**2/(arg['umax1']*arg['m2Z']*arg['kZ1']**2)
    cpo4 = arg['kP1']/arg['umax1']*(arg['mP1']-arg['m1Z']*arg['gmax1']/(arg['m2Z']*arg['kZ1']))
                                    
    a = arg['epsilon2Z']*(1-arg['gamma'])*arg['gmax1']/arg['kZ1']*az-apo4*arg['umax1']/arg['kP1']
    b = -arg['epsilon2Z']*(1-arg['gamma'])*arg['gmax1']/arg['kZ1']*cz+arg['mP1']-cpo4*arg['umax1']/arg['kP1']+arg['epsilon1Z']*arg['m1Z']*az
    c = p-arg['epsilon1Z']*arg['m1Z']*cz
    delta = b**2-4*a*c 
    x1 = (-b-math.sqrt(delta))/(2*a)
    x2 =(-b+math.sqrt(delta))/(2*a)
    
    P1_barre = x1
    P2_barre = 0
    Z_barre = az*P1_barre-cz
    PO4_barre = apo4*P1_barre+cpo4

if Phyto == 'P1null':
    
    p = 0.1
    Psupply = [p] * len(time)
    
    [P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10_simplify(Psupply, time, P1_ini = 0)
    
    az = arg['gamma']*arg['gmax2']/(arg['m2Z']*arg['kZ2'])
    cz = arg['m1Z']/arg['m2Z']
    apo4 = arg['kP2']*arg['gamma']*arg['gmax2']**2/(arg['umax2']*arg['m2Z']*arg['kZ2']**2)
    cpo4 = arg['kP2']/arg['umax2']*(arg['mP2']-arg['m1Z']*arg['gmax2']/(arg['m2Z']*arg['kZ2']))                                
    a = arg['epsilon2Z']*(1-arg['gamma'])*arg['gmax2']/arg['kZ2']*az-apo4*arg['umax2']/arg['kP2']
    b = -arg['epsilon2Z']*(1-arg['gamma'])*arg['gmax2']/arg['kZ2']*cz+arg['mP2']-cpo4*arg['umax2']/arg['kP2']+arg['epsilon1Z']*arg['m1Z']*az
    c = p-arg['epsilon1Z']*arg['m1Z']*cz  
    delta = b**2-4*a*c 
    x1 = (-b-math.sqrt(delta))/(2*a)
    x2 =(-b+math.sqrt(delta))/(2*a)
       
    P1_barre = 0
    P2_barre = x1
    Z_barre = az*P2_barre-cz
    PO4_barre = apo4*P2_barre+cpo4
         
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