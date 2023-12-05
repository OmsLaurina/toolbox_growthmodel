"""
Set up function
"""

import numpy as np

def set_up_steadystate():
    
    dt = 0.1 # time step
    end_time = 2000 # end of simulation time
    time = np.arange(0, end_time, dt) # simulation time 
    nb_time = len(time) # number if time step
    
    Psupply_cst = 0.1 # Constant Psupply 
    Psupply = [Psupply_cst] * nb_time # Psupply list
    
    Psupply_1 = 1 
    Psupply_senstest = [Psupply_1] * nb_time
    Psupply_senstest = np.array(Psupply_senstest)
    
    n = 100 # nb of tested values
    min_param = 0.01 # min 
    max_param = 0.1  # max
    l_param = np.linspace(min_param,max_param,n) 

    return dt, end_time, time, Psupply, Psupply_cst, Psupply_senstest, min_param, max_param,l_param

def set_up_dynamicalstate():
    
    dt = 0.1 # time step
    end_time_flux = 90 # end of simulation time
    time_flux = np.arange(0,end_time_flux,dt)
    nb_time_flux=len(time_flux)
    
    Psupply_ini = 0.05 # psupply at t=0
    b = 0.08 # amplitude
    nbpulse = 1 # number of pulse

    return dt, Psupply_ini, b, nbpulse, end_time_flux, time_flux, nb_time_flux

def set_up_diagrbifurc():
    
    dt = 0.1
    end_time = 2000
    time = np.arange(0,end_time,dt)
    Psupply_moy = 1
    Psupply_arr = np.array([Psupply_moy]*len(time))
    
    n = 100  
    min_param = 0.01
    max_param = 0.4
    l_param = np.linspace(min_param,max_param,n)
    
    return dt, end_time, time, Psupply_moy, Psupply_arr, n, l_param
    
