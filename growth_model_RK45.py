import numpy as np

# Set up

# Define the function of growth model

def growth_model(t, state, Psupply):
      
    P1, P2, Z, PO4 = state
    
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
        }
    
# Growth rates (monod function)
    u1 = (PO4 / (arg['kP1'] + PO4)) * arg['umax1']
    u2 = (PO4 / (arg['kP2'] + PO4)) * arg['umax2']

    # Functional response = grazing rates (holling type II)
    g1 = (P1 / (arg['kZ1'] + P1 + P2)) * arg['gmax1']
    g2 = (P2 / (arg['kZ2'] + P1 + P2)) * arg['gmax2']

    # Fluxes
    PP1 = u1 * P1
    PP2 = u2 * P2
    G1 = g1 * Z
    G2 = g2 * Z
    exc = (1 - arg['gamma']) * g1 * Z + (1 - arg['gamma']) * g2 * Z
    d_P1 = arg['mP1'] * P1
    d_P2 = arg['mP2'] * P2
    d1_Z = arg['m1Z'] * Z
    d2_Z = arg['m2Z'] * Z ** 2
    rec_P1 = arg['epsilonP'] * d_P1
    rec_P2 = arg['epsilonP'] * d_P2
    rec_Z = arg['epsilon1Z'] * d1_Z
    rec_exc = arg['epsilon2Z'] * exc

    # P1 takes from the nutrient pool first
    max_available_PO4 = PO4 +  rec_P1  + rec_P2  + rec_Z  + rec_exc 

    # if PP1 > max_available_PO4:
    #     PP1 = max_available_PO4

    # max_available_PO4 = max_available_PO4 - PP1 
    # if PP2 > max_available_PO4:
    #     PP2 = max_available_PO4

    # Derivatives
    dP1_dt = PP1 - G1 - d_P1
    dP2_dt = PP2 - G2 - d_P2
    dZ_dt = arg['gamma'] * (G1 + G2) - d1_Z - d2_Z
    dPO4_dt =  Psupply + rec_P1 + rec_P2 + rec_Z  + rec_exc  - PP1  - PP2 
    
    return [dP1_dt, dP2_dt, dZ_dt, dPO4_dt] 