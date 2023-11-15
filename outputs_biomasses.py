import numpy as np
import pandas as pd

# Load data of abundances and FWS
index_hipp = slice(408, 423) # index for the transect showed by Tzortzis et al 2021: 408, 423, and for the full hippodrome: 408, 511
n = index_hipp.stop - index_hipp.start

data_ab = np.loadtxt('../BioSWOTmed_data/data_CYTO_NEW.txt', delimiter=';')
data_ab = data_ab[index_hipp, :]

# For the full hippodrome 
# ligne_a_supprimer = [16, 17]
# data_ab = np.delete(data_ab, ligne_a_supprimer, axis=0)

lat_ab = data_ab[:, 23]
lon_ab = data_ab[:, 22]

mask = lat_ab > 38.5

### Extract abundances from files and specific parameters for the chosen group

chosen_group = 'MICRO'

if chosen_group == 'SYNE':
    ab = data_ab[:, 6]
    data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL synechococcus .txt', delimiter=';') 
    
    # log-log regression coefficients from Menden-Deurer 2000 Table.4: 
    a = 0.21
    b = 0.939

if chosen_group == 'PICO1':
    ab = data_ab[:, 8]
    data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL red_pico_euk .txt', delimiter=';')
    
    # log-log regression coefficients from Menden-Deurer 2000 Table.4: 
    a = 0.21
    b = 0.939

if chosen_group == 'PICO2':
    ab = data_ab[:, 9]
    data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL red_pico_2 .txt', delimiter=';')
    
    # log-log regression coefficients from Menden-Deurer 2000 Table.4: 
    a = 0.21
    b = 0.939

if chosen_group == 'PICO3':
    ab = data_ab[:, 10]
    data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL red-pico-euk-3 .txt', delimiter=';')
    data_FWS = data_FWS.replace(',', '.', regex=True)
    
    # To convert text to number (eg: 1E3 --> 10³)
    for colonne in data_FWS.columns:
        try:
            data_FWS[colonne] = pd.to_numeric(data_FWS[colonne], errors='coerce')
        except ValueError:
            pass
        
    # log-log regression coefficients from Menden-Deurer 2000 Table.4: 
    a = 0.21
    b = 0.939

if chosen_group == 'PICOHFLR':
    ab = data_ab[:, 14]
    data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL pico_HFLR .txt', delimiter=';')
    data_FWS = data_FWS.replace(',', '.', regex=True)
    
    # To convert text to number (eg: 1E3 --> 10³)
    for colonne in data_FWS.columns:
        try:
            data_FWS[colonne] = pd.to_numeric(data_FWS[colonne], errors='coerce')
        except ValueError:
            pass
    
    # log-log regression coefficients from Menden-Deurer 2000 Table.4: 
    a = 0.21
    b = 0.939

if chosen_group == 'NANOsws':
    ab = data_ab[:, 11]
    data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL red-nano-euk-sws .txt', delimiter=';')
    
    # log-log regression coefficients from Menden-Deurer 2000 Table.4: 
    a = 0.26
    b = 0.86

if chosen_group == 'NANOred':
    ab = data_ab[:, 12]
    data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL red_nano_euk .txt', delimiter=';')
    
    # log-log regression coefficients from Menden-Deurer 2000 Table.4: 
    a = 0.26
    b = 0.86

if chosen_group == 'MICRO':
    ab = data_ab[:, 13]
    data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL microphyto .txt', delimiter=';')
    
    # log-log regression coefficients from Menden-Deurer 2000 Table.4: 
    a = 0.26
    b = 0.86
    
if chosen_group == 'CRYPTO':
    ab = data_ab[:, 7]
    data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL cryptophyte .txt', delimiter=';')
    
    # log-log regression coefficients from Menden-Deurer 2000 Table.4: 
    a = 0.26
    b = 0.86

# Transform table
FWS = data_FWS.iloc[:, 0:n].values

### Convert abundances to biomasses

## Calculate bioV [um³]

# log-log regression coefficient bioV(FWS) (cf. Tzortzis et al 2023)
beta_0 = -5.8702
beta_1 = 0.9228

bioV_um3 = FWS**beta_1 * np.exp(beta_0)

# mean bioV of the chosen group
bioV_um3 = np.nanmean(bioV_um3, axis=0) 
bioV_north = bioV_um3[mask]
bioV_south = bioV_um3[~mask]

## Calculate Qc [pgC/cell]

Qc = a*bioV_um3**b

# [fgC/cell]
Qc_north1 = (a*bioV_north**b)*1e3
Qc_south1 = (a*bioV_south**b)*1e3

## Calculate biomass biom [mmolC/m³]

Qc = Qc * 1e-12  # pg --> g
Qc = Qc / 12.011  # g --> molC
Qc = Qc * 1e3  # molC --> mmolC

ab = ab * 1e6  # cm3 --> m3
ab_north = ab[mask]
ab_south = ab[~mask]

Qc_north = Qc.T[mask]
Qc_south = Qc.T[~mask]

biom_north = ab_north * np.nanmean(Qc_north)
biom_south = ab_south * np.nanmean(Qc_south)
biom = ab*np.nanmean(Qc)

### Save biomasses and corresponding lon/lat

np.savetxt(f'../outputs/{chosen_group}.txt', biom)
np.savetxt('../outputs/lon.txt', lon_ab)
np.savetxt('../outputs/lat.txt', lat_ab)

### Create table with mean values of each variables

# ["Abundance [cell/cm3]", "Biovolume [um3/cell]", "Qc [fgC/cell]", "Biomass [mmolC/m3]"]

T_north = [np.mean(ab_north)/1e6, np.nanmean(bioV_north), np.nanmean(Qc_north1), np.mean(biom_north)]
T_south = [np.mean(ab_south)/1e6, np.nanmean(bioV_south), np.nanmean(Qc_south1), np.mean(biom_south)]

T_north = [
    format(np.mean(ab_north)/1e6, ".4e"),
    round(np.nanmean(bioV_north), 4),
    format(np.nanmean(Qc_north1), ".4e"),
    round(np.nanmean(biom_north), 4)
]

T_south = [
    format(np.mean(ab_south)/1e6, ".4e"),
    round(np.nanmean(bioV_south), 4),
    format(np.nanmean(Qc_south1), ".4e"),
    round(np.mean(biom_south), 4)
]

merged_array = np.array([T_north, T_south])
merged_array = merged_array.T

rows = ["Abundance [cell/cm3]", "Biovolume [um3/cell]", "Qc [fgC/cell]", "Biomass [mmolC/m3]"]
columns = ["North water", "South water"]
df = pd.DataFrame(merged_array, columns=columns, index=rows)

df.to_csv(f'../outputs/table_{chosen_group}.csv', sep='\t')