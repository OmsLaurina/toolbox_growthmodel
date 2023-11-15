""" PLOT CHLOROPHYLL DURING PROTEVSWOT 2018 - MEDSEA_MULTIYEAR_BGC_006_008 """

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
from matplotlib.colors import Normalize
import pandas as pd

plt.close('all')

file_path = "med-ogs-pft-rean-d_1698681358118.nc"
ds = xr.open_dataset(file_path)

chl = ds['chl']
lat = ds['latitude']
lon = ds['longitude']

output_directory = "../figures"

# Sélectionner la plage de dates du 8 mai au 15 mai
start_date = '2018-04-30T12:00:00'
end_date = '2018-05-15T12:00:00'

chlorophylle_maps = []

for single_date in pd.date_range(start=start_date, end=end_date, freq='D'):
    chl_in_date = chl.sel(time=single_date)
    
    plt.figure(figsize=(6, 4))
    im = chl_in_date.plot(x='longitude', y='latitude', cmap='viridis')
    norm = Normalize(vmin=0.05, vmax=0.1)
    im.set_norm(norm)
    plt.ylim([37.5, 40.5])
    plt.xlim([2, 6])
    plt.title(f'Chlorophylle du {single_date}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    output_file = os.path.join(output_directory, f'chlorophylle_map_{single_date}.png')
    plt.savefig(output_file)
    plt.close()

plt.show()


    

