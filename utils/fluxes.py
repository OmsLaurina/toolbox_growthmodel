import numpy as npy
import numpy as np
import random

def periodicflux(Psupply_moy,time_flux,**kwargs): #Grover 1990, Mayersohn 2022
    global Psupply_periodic
    nb_time=len(time_flux)
    Psupply = np.zeros(nb_time)
    defaultKwargs = {
        'b':0.01, # Amplitude
        'T':30,   # Periode
        }
    arg = { **defaultKwargs, **kwargs}
    for k in range(0,nb_time):
        Psupply[k] = Psupply_moy+arg['b']*np.sin(2*np.pi*k/arg['T'])
    return Psupply,arg

def pulsedfluxGrover1990(Psupply_moy,time_flux,**kwargs): #Grover 1990
    global Psupply_pulsed
    nb_time=len(time_flux)
    Psupply = np.zeros(nb_time)
    defaultKwargs = {
        'D':0.25        }
    arg = { **defaultKwargs, **kwargs}
    min_T = 5
    max_T = 30
    l_T = np.linspace(min_T, max_T, nb_time)
    b_list = []
    for i in range(0,nb_time):
        b = Psupply_moy * arg['D'] * l_T[i]/ (1 - np.exp(arg['D']) * l_T[i]) - Psupply_moy
        b_list.append(b)
    for k in range(0,nb_time):
        Psupply[k] = Psupply_moy+b_list[k]*np.sin(2*np.pi*k/l_T[k])
    return Psupply,arg


def pulsedflux_randomAmplitude(Psupply_moy,time_flux,**kwargs): #Random amplitude
    global Psupply_pulsed
    nb_time=len(time_flux)
    Psupply = np.zeros(nb_time)
    defaultKwargs = {
        'T':5,   # Periode
        }
    arg = { **defaultKwargs, **kwargs}
    b = [random.uniform(0.01, 0.1) for _ in range(nb_time)]
    for k in range(0,nb_time):
        Psupply[k] = Psupply_moy+b[k]*np.sin(2*np.pi*k/arg['T'])
    return Psupply,arg


def pulsedflux_stepfunction(time_flux, **kwargs):
    nb_time = len(time_flux)
    Psupply = np.zeros(nb_time)
    default_kwargs = {
        'A': -0.1,
        'C': 0.01
    }
    arg = {**default_kwargs, **kwargs}
    pulsations = kwargs.get('pulsations', [])
    for pulsation in pulsations:
        t1 = pulsation['t1']
        t2 = pulsation['t2']
        Psupply += arg['A'] * (np.heaviside(time_flux - t1, 1) - np.heaviside(time_flux - t2, 1))
    Psupply += arg['C']
    return Psupply,arg