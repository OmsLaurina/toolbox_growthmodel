#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 17:38:17 2023

@author: loms
"""

import numpy as np
import matplotlib.pyplot as plt

biom_Syne = np.loadtxt('../outputs/Syne.txt')
biom_Micro = np.loadtxt('../outputs/MICRO.txt')
biom_PICO1 = np.loadtxt('../outputs/PICO1.txt')
biom_PICO2 = np.loadtxt('../outputs/PICO2.txt')
biom_PICO3 = np.loadtxt('../outputs/PICO3.txt')
biom_NANOred = np.loadtxt('../outputs/NANOred.txt')
biom_NANOsws = np.loadtxt('../outputs/NANOsws.txt')
biom_PICOHFLR = np.loadtxt('../outputs/PICOHFLR.txt')
lat = np.loadtxt('../outputs/lat.txt')

mask = lat > 38.5

tot = biom_Syne+biom_Micro+biom_NANOred+biom_NANOsws+biom_PICO1+biom_PICO3+biom_PICO2+biom_PICOHFLR
pico = biom_Syne+biom_PICO1+biom_PICO2+biom_PICO3+biom_PICOHFLR
nano = biom_NANOred+biom_NANOsws
micro = biom_Micro 

np.savetxt('../outputs/tot.txt', tot)
np.savetxt('../outputs/picob.txt', pico)
np.savetxt('../outputs/nanob.txt', nano)
np.savetxt('../outputs/microb.txt', micro)

# prop_pico = (pico/tot)*100
prop_pico = pico
prop_pico_low = prop_pico[~mask]
prop_pico_high = prop_pico[mask]
np.savetxt('../outputs/prop_pico.txt', prop_pico)

# prop_nano = (nano/tot)*100
prop_nano = nano
prop_nano_low = prop_nano[~mask]
prop_nano_high = prop_nano[mask]
np.savetxt('../outputs/prop_nano.txt', prop_nano)

# prop_micro = (micro/tot)*100
prop_micro = micro
prop_micro_low = prop_micro[~mask]
prop_micro_high = prop_micro[mask]
np.savetxt('../outputs/prop_micro.txt', prop_micro)

f_biomass_syne = np.mean(prop_pico_low)/np.mean(prop_pico_high)
f_biomass_nano = np.mean(prop_nano_low)/np.mean(prop_nano_high)
f_biomass_micro = np.mean(prop_micro_low)/np.mean(prop_micro_high)
factor = np.column_stack((f_biomass_syne,f_biomass_nano,f_biomass_micro))
np.savetxt('../outputs/f_biomass.txt', factor, header = 'f_syne f_nano f_micro')


R = pico/np.mean(tot)
Rsud = R[mask]
