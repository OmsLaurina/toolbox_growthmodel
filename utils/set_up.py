"""
Set up function
"""

import numpy as np

def set_up_steadystate():
    dt = 0.1
    end_time = 2000
    time = np.arange(0, end_time, dt)
    Psupply_cst = 0.05
    nb_time = len(time)
    Psupply = [Psupply_cst] * nb_time
    
    # Psupply = 1 for the sensitivity test
    Psupply_1 = 1
    Psupply_senstest = [Psupply_1] * nb_time
    Psupply_senstest = np.array(Psupply_senstest)
    
    # list of different value of psupply for the sensitivity test
    n = 100 # nb of tested values
    min_param = 0.01 # min 
    max_param = 0.1  # max
    l_param = np.linspace(min_param,max_param,n)

    return dt, end_time, time, Psupply, Psupply_cst, Psupply_senstest, l_param
    
