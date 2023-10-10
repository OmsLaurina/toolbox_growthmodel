"""

Record the in situ biomasses and calculate the fbioM ratio. 

"""

import numpy as np
import matplotlib.pyplot as plt

# Load data from the BioSWOTmed_data directory
biom_Syne = np.loadtxt('../BioSWOTmed_data/SYNE.txt')
biom_Micro = np.loadtxt('../BioSWOTmed_data/MICRO.txt')
biom_PICO1 = np.loadtxt('../BioSWOTmed_data/PICO1.txt')
biom_PICO2 = np.loadtxt('../BioSWOTmed_data/PICO2.txt')
biom_PICO3 = np.loadtxt('../BioSWOTmed_data/PICO3.txt')
biom_NANOred = np.loadtxt('../BioSWOTmed_data/NANOred.txt')
biom_NANOsws = np.loadtxt('../BioSWOTmed_data/NANOsws.txt')
biom_PICOHFLR = np.loadtxt('../BioSWOTmed_data/PICOHFLR.txt')
lat = np.loadtxt('../BioSWOTmed_data/lat.txt')

# Calculate and record the regrouped biomasses
tot = biom_Syne+biom_Micro+biom_NANOred+biom_NANOsws+biom_PICO1+biom_PICO3+biom_PICO2+biom_PICOHFLR
pico = biom_Syne+biom_PICO1+biom_PICO2+biom_PICO3+biom_PICOHFLR
nano = biom_NANOred+biom_NANOsws
micro = biom_Micro 

np.savetxt('../outputs/tot.txt', tot)
np.savetxt('../outputs/pico.txt', pico)
np.savetxt('../outputs/nano.txt', nano)
np.savetxt('../outputs/micro.txt', micro)

# Create a mask at the front position to extract south and north biomass
mask = lat > 38.5

# pico = (pico/tot)*100
pico_low = pico[~mask]
pico_high = pico[mask]
# nano = (nano/tot)*100
nano_low = nano[~mask]
nano_high = nano[mask]
# micro = (micro/tot)*100
micro_low = micro[~mask]
micro_high = micro[mask]

# Calculate the fbioM ratio
f_biomass_syne = np.mean(pico_low)/np.mean(pico_high)
f_biomass_nano = np.mean(nano_low)/np.mean(nano_high)
f_biomass_micro = np.mean(micro_low)/np.mean(micro_high)
factor = np.column_stack((f_biomass_syne,f_biomass_nano,f_biomass_micro))
np.savetxt('../outputs/f_biomass.txt', factor, header = 'f_syne f_nano f_micro')