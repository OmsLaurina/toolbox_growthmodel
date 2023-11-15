import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np
import xarray as xr
from matplotlib.colors import Normalize
import scipy.io

plt.close('all')


# Load data of abundances and FWS
index_hipp = slice(408, 511)  # 408 to 511 pour l'hipp complet, 408 to 420 pour le transect de Roxane

data_ab = np.loadtxt('../BioSWOTmed_data/data_CYTO_NEW.txt', delimiter=';')
data_ab = data_ab[index_hipp, :] 
# ligne_a_supprimer = [16, 17]
# data_ab = np.delete(data_ab, ligne_a_supprimer, axis=0)

lat_ab = data_ab[:, 23]
lon_ab = data_ab[:, 22]

# Liste des groupes
group_names = ['Syne', 'PICO1', 'PICO2', 'PICO3', 'PICOHFLR', 'NANOsws', 'NANOred', 'MICRO', 'All']

#fig, axes = plt.subplots(2, 4, figsize=(15, 15))

# file_path = "med-ogs-pft-rean-d_1698681358118.nc"
# ds = xr.open_dataset(file_path)

# chl = ds['chl']
# lat = ds['latitude']
# lon = ds['longitude']
# time = ds['time']

# start_date = '2018-04-30T12:00:00'
# end_date = '2018-05-18T12:00:00'
# chl_in_date_range = chl.sel(time=slice(start_date, end_date))
# chl_mean = chl_in_date_range.mean(dim='time')

chl_data = scipy.io.loadmat('../BioSWOTmed_data/nrt_color_swmed_20180510_20180512.mat')
lon = chl_data['lon']
lat = chl_data['lat']
chl = chl_data['chl']
lon, lat = np.meshgrid(lon, lat)

for i, chosen_group in enumerate(group_names):
    plt.figure(figsize=(6, 4))

    if chosen_group == 'Syne':
        ab = data_ab[:, 6]
        data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL synechococcus .txt', delimiter=';') 
        a = 0.26
        b = 0.86

    if chosen_group == 'PICO1':
        ab = data_ab[:, 8]
        data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL red_pico_euk .txt', delimiter=';')
        a = 0.26
        b = 0.86

    if chosen_group == 'PICO2':
        ab = data_ab[:, 9]
        data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL red_pico_2 .txt', delimiter=';')
        a = 0.26
        b = 0.86

    if chosen_group == 'PICO3':
        ab = data_ab[:, 10]
        data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL red-pico-euk-3 .txt', delimiter=';')
        data_FWS = data_FWS.replace(',', '.', regex=True)
        for colonne in data_FWS.columns:
            try:
                data_FWS[colonne] = pd.to_numeric(data_FWS[colonne], errors='coerce')
            except ValueError:
                pass
        a = 0.26
        b = 0.86

    if chosen_group == 'PICOHFLR':
        ab = data_ab[:, 14]
        data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL pico_HFLR .txt', delimiter=';')
        data_FWS = data_FWS.replace(',', '.', regex=True)
        for colonne in data_FWS.columns:
            try:
                data_FWS[colonne] = pd.to_numeric(data_FWS[colonne], errors='coerce')
            except ValueError:
                pass
        a = 0.26
        b = 0.86

    if chosen_group == 'NANOsws':
        ab = data_ab[:, 11]
        data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL red-nano-euk-sws .txt', delimiter=';')
        a = 0.433
        b = 0.863

    if chosen_group == 'NANOred':
        ab = data_ab[:, 12]
        data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL red_nano_euk .txt', delimiter=';')
        a = 0.433
        b = 0.863

    if chosen_group == 'MICRO':
        ab = data_ab[:, 13]
        data_FWS = pd.read_csv('../BioSWOTmed_data/FWS_ALL microphyto .txt', delimiter=';')
        a = 0.433
        b = 0.863
        
    if chosen_group == 'All':
        ab = data_ab[:, 6]+data_ab[:, 8]+data_ab[:, 9]+data_ab[:, 10]+data_ab[:, 11]+data_ab[:, 12]+data_ab[:, 13]+data_ab[:, 14]
        
    # im = chl_mean.plot(x='longitude', y='latitude', cmap='plasma',cbar_kwargs={"shrink": 0.5})
    im = plt.imshow(chl, extent=(lon.min(), lon.max(), lat.min(), lat.max()), origin='lower', cmap='viridis')
    plt.colorbar(im, label=r'chlorophylle $[mg.m^{-3}]$',shrink=0.5)
    
    # Définir une échelle de couleur personnalisée
    norm = Normalize(vmin=0.1, vmax=0.2)
    im.set_norm(norm)
    # plt.colorbar(label='Chlorophylle')
    # norm = Normalize(vmin=0.05, vmax=0.1)
    # im.set_norm(norm)

    scatter = plt.scatter(lon_ab, lat_ab, c=ab, cmap=cm.plasma, s=20)
    cbar2 = plt.colorbar(scatter, shrink=0.5)
    cbar2.set_label(r'abundance $[cell.cm^{-3}]$', rotation=90)

    plt.title(chosen_group)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.ylim([37.5, 39.5]) 
    plt.xlim([2, 4]) 

    # # Sélectionner le subplot actuel
    # row = i // 4
    # col = i % 4
    # ax = axes[row, col]
    
    # im = chl_on_target_date.plot(x='longitude', y='latitude', cmap='plasma', ax=ax)
    
    # cbar = plt.colorbar(im, ax=ax, shrink=0.5)
    # cbar.set_label('chl', rotation=90)
    # norm = Normalize(vmin=0.05, vmax=0.1)
    # im.set_norm(norm)

    # scatter = ax.scatter(lon_ab, lat_ab, c=ab, cmap=cm.viridis, s=40)
    
    # cbar2 = plt.colorbar(scatter, ax=ax, shrink=0.5)
    # cbar2.set_label('ab', rotation=90)
    
    # ax.set_title(chosen_group)
    # ax.set_xlabel('Longitude')
    # ax.set_ylabel('Latitude')
    # ax.set_ylim([37.5, 40.5]) 
    # ax.set_title(chosen_group)

    plt.tight_layout()
    plt.show()
    plt.savefig(f'../figures/figures_abundancesfull_{chosen_group}.pdf', format='pdf')