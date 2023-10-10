import numpy as np
from growth_model_2P1Z_v10 import growth_model_2P1Z_v10
from fluxes import pulsedflux_stepfunction
import matplotlib.pyplot as plt
import numpy as np
import scienceplots
plt.close('all')

# Timestep
dt = 0.1
end_time = 2000
time = np.arange(1, end_time, dt)

# External input of nutrients
Psupply_ini = 0.05
moy_flux = Psupply_ini
Psupply_cst = np.full(len(time), Psupply_ini)
b = 0.08
nbpulse = 1

end_time_flux = 2000
time_flux = np.arange(1, end_time_flux, dt)

fluxes = 'pulsedflux_stepfunction'

if nbpulse == 3:
    [Psupply, arg] = pulsedflux_stepfunction(time_flux, A=b, C=moy_flux, pulsations=[{'t1': 5, 't2': 10}, {'t1': 20, 't2': 25}, {'t1': 35, 't2': 40}])
elif nbpulse == 2:
    [Psupply, arg] = pulsedflux_stepfunction(time_flux, A=b, C=moy_flux, pulsations=[{'t1': 5, 't2': 10}, {'t1': 20, 't2': 25}])
elif nbpulse == 1:
    [Psupply, arg] = pulsedflux_stepfunction(time_flux, A=b, C=moy_flux, pulsations=[{'t1': 5, 't2': 10}])

# Get the solution at steady state equilibrium
[P1, P2, Z, PO4, arg] = growth_model_2P1Z_v10(Psupply_cst, time_flux, dt)
P1_eq, P2_eq, Z_eq, PO4_eq = P1[-1], P2[-1], Z[-1], PO4[-1]

#10% autour de la valeur d'equilibre
tx = 0.1

P1_eq_var = P1_eq*(1+tx)
P2_eq_var = P2_eq*(1+tx)
Z_eq_var = Z_eq*(1+tx)
PO4_eq_var = PO4_eq*(1+tx)
V1 = np.sqrt((P1_eq_var-P1_eq)**2+(P2_eq_var-P1_eq)**2+(Z_eq_var-Z_eq)**2+(PO4_eq_var-PO4_eq)**2)
C1 = (P1_eq_var-P1_eq)**2+(P2_eq_var-P1_eq)**2+(Z_eq_var-Z_eq)**2+(PO4_eq_var-PO4_eq)**2

# New solutions
P1, P2, Z, PO4, _ = growth_model_2P1Z_v10(Psupply, time_flux, dt, P1_ini=P1_eq, P2_ini=P2_eq, Z_ini=Z_eq, PO4_ini=PO4_eq) 
  
bo = True
# ### Method 1 : calcul de la sphère autour de l'equilibre
# epsilon_carre = 10**(-3.5)
epsilon_carre = []
for i, t in enumerate(time_flux):
    if t>10:
        C = (P1[i]-P1_eq)**2+(P2[i]-P2_eq)**2+(Z[i]-Z_eq)**2+(PO4[i]-PO4_eq)**2
        epsilon_carre = abs(C1-C)
        if C<epsilon_carre and bo:
            bo = False
            t_res = t
            print('resilience time:', t)

### Method 2 : calcul de la norme du vecteur  
# epsilon = []
# # epsilon = 10**(-6)       
# for i, t in enumerate(time_flux[0:len(time_flux)-1]):
#     if t>10:
#         V = np.sqrt((P1[i+1]-P1[i])**2+(P2[i+1]-P2[i])**2+(Z[i+1]-Z[i])**2+(PO4[i+1]-PO4[i])**2)
#         espilon = abs(V1-V)
#         if V<epsilon and bo:
#             bo = False
#             t_res = t
#             print('resilience time:', t)         
            
plt.figure(1)
plt.rc('font', size=7)
plt.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time_flux, P2, label=r'$P_2$', color="green")
plt.plot(time_flux, Z, label=r'$Z$', color="aqua")
plt.plot(time_flux, PO4, label=r'$PO4$', color="magenta")
plt.plot(time_flux, Psupply, label='Psupply',color="gray")
plt.axvline(x=t_res, color='red', linestyle='--')
plt.legend()