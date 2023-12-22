"""

Plot the in situ biomasses. Create the figure 2 of the paper and the left panel of the figure 9.

"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker
import scienceplots
from scipy.interpolate import griddata
import xarray as xr
import scipy.io
from scipy.interpolate import interp2d
from sklearn.preprocessing import MinMaxScaler

plt.style.use(['science','no-latex'])
plt.close('all')

# Load data from outputs_insitu.py code
pico = np.loadtxt('../outputs/pico.txt')
nano = np.loadtxt('../outputs/nano.txt')
micro = np.loadtxt('../outputs/micro.txt')
nanosws = np.loadtxt('../outputs/NANOsws.txt')
nanored = np.loadtxt('../outputs/NANOred.txt')
crypto = np.loadtxt('../outputs/CRYPTO.txt')
tot = np.loadtxt('../outputs/tot.txt')

lat = np.loadtxt('../outputs/lat.txt')
lon = np.loadtxt('../outputs/lon.txt')

f_BioM = np.loadtxt('../outputs/f_biomass.txt')
f_biomass_syne,f_biomass_nano,f_biomass_micro = f_BioM.T

# Figure 1 (figure 2 of the paper): Biomasses distribution
fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(1, 3, 1) # PICO
ax2 = fig.add_subplot(1, 3, 2) # MICRO
ax3 = fig.add_subplot(1, 3, 3) # NANO
# ax4 = fig.add_subplot(1, 4, 4) # TOT

normal_size = 20

# scatter1 = ax1.scatter(lon, lat, c=pico, cmap='viridis')
scatter1 = ax1.scatter(lon, lat, c=nanored, s=[normal_size if i > 13 else normal_size*3.5 for i in range(len(lon))], cmap='viridis')
ax1.set_title(r'$PICO$' + '\n$f_{BioM}=' + str(round(f_biomass_syne, 1)) + '$')
ax1.set_xlabel('Longitude [°E]')
ax1.set_ylabel('Latitude [°N]')
ax1.set_aspect('equal')
cbar1 = plt.colorbar(scatter1, ax=ax1, label='Biomass [mmolC m$^{-3}$]')
ax1.axhline(y=38.5, color='lightcoral', linestyle='-',linewidth=1.5)
tick_locator = ticker.MaxNLocator(nbins=10)
cbar1.locator = tick_locator
cbar1.update_ticks()

# scatter2 = ax2.scatter(lon, lat, c=nano, cmap='viridis')
scatter2 = ax2.scatter(lon, lat, c=nano, s=[normal_size if i > 13 else normal_size*3.5 for i in range(len(lon))], cmap='viridis')
ax2.set_title(r'$NANO$' + '\n$f_{BioM}=' + str(round(f_biomass_nano, 1)) + '$')
ax2.set_xlabel('Longitude [°E]')
ax2.set_ylabel('Latitude [°N]')
ax2.set_aspect('equal')
cbar2 = plt.colorbar(scatter2, ax=ax2, label='Biomass [mmolC m$^{-3}$]')
ax2.axhline(y=38.5, color='lightcoral', linestyle='-',linewidth=1.5)
tick_locator = ticker.MaxNLocator(nbins=10)
cbar2.locator = tick_locator
cbar2.update_ticks()

# scatter3 = ax3.scatter(lon, lat, c=micro, cmap='viridis')
scatter3 = ax3.scatter(lon, lat, c=micro, s=[normal_size if i > 13 else normal_size*3.5 for i in range(len(lon))], cmap='viridis')
ax3.set_title(r'$MICRO$' + '\n$f_{BioM}=' + str(round(f_biomass_micro, 1)) + '$')
ax3.set_xlabel('Longitude [°E]')
ax3.set_ylabel('Latitude [°N]')
ax3.set_aspect('equal')
cbar3 = plt.colorbar(scatter3, ax=ax3, label='Biomass [mmolC m$^{-3}$]')
ax3.axhline(y=38.5, color='lightcoral', linestyle='-',linewidth=1.5)
tick_locator = ticker.MaxNLocator(nbins=10)
cbar3.locator = tick_locator
cbar3.update_ticks()

# scatter4 = ax4.scatter(lon, lat, c=tot, cmap='viridis')
# ax4.set_title(r'$TOTAL$')
# ax4.set_xlabel('Longitude [°E]')
# ax4.set_ylabel('Latitude [°N]')
# ax4.set_aspect('equal')
# cbar4 = plt.colorbar(scatter4, ax=ax4, label='Biomass [mmolC m$^{-3}$]')
# ax4.axhline(y=38.5, color='lightcoral', linestyle='-',linewidth=1.5)
# tick_locator = ticker.MaxNLocator(nbins=10)
# cbar4.locator = tick_locator
# cbar4.update_ticks()

plt.tight_layout()
plt.show()
plt.savefig('../figures/insitu_biomasses.pdf', format='pdf')

# Calculate the polynomial regression of the in situ R ratio (pico/mean(to))
degree = 3
coefficients = np.polyfit(lat, (pico+nanored)/(micro+pico+nano), degree)
polynomial = np.poly1d(coefficients)
x_fit = np.linspace(min(lat), max(lat), 1000)
y_fit = polynomial(x_fit)

# Figure 2 (figure 9 of the paper): R ratio as a function of the latitude
plt.figure(2)
plt.scatter(lat, (pico+nanored)/(micro+pico+nano), c=(pico+micro+nano), cmap='Wistia')
plt.plot(x_fit, y_fit, 'k-', label=f'Polynomial Fit (Degree {degree})')
plt.axvline(x=38.5, color='red', linestyle='--',linewidth=2)
plt.grid()
plt.colorbar(label=r'${TotBioM}$ [mmolC m$^{-3}$]')
plt.xlabel('Latitude [°N]', fontsize=10)
plt.ylabel(r'$R_{in situ}$ []', fontsize=10)
plt.savefig('../figures/ratioR_latitude.pdf', format='pdf')

# Same than the previous one but with scatter size depending on chlorophyll

chl_data = scipy.io.loadmat('../BioSWOTmed_data/chl_CLS/nrt_color_swmed_20180510_20180512.mat')
lon_chl = chl_data['lon'].T
lat_chl = chl_data['lat'].T
chl = chl_data['chl']

lat_biomasses = np.loadtxt('../outputs/lat.txt')
lon_biomasses = np.loadtxt('../outputs/lon.txt')

lon_chl_adjusted, lat_chl_adjusted = np.meshgrid(lon_chl.flatten(), lat_chl.flatten())
points = np.column_stack((lon_chl_adjusted.flatten(), lat_chl_adjusted.flatten()))
values = chl.flatten()
chl_interpolated = griddata(points, values, (lon_biomasses, lat_biomasses), method='linear', fill_value=np.nan)

scaler = MinMaxScaler()
chl_at_lat_normalized = scaler.fit_transform(chl_interpolated.reshape(-1, 1)).flatten()
size = chl_at_lat_normalized * 100

plt.figure(3)
plt.scatter(lat, (pico+nanored)/(micro+pico+nano), c=(pico+micro+nano), cmap='Wistia', s=size)
plt.plot(x_fit, y_fit, 'k-', label=f'Polynomial Fit (Degree {degree})')
plt.axvline(x=38.5, color='red', linestyle='--', linewidth=2)
plt.grid()
plt.colorbar(label=r'${TotBioM}$ [mmolC m$^{-3}$]')
plt.xlabel('Latitude [°N]', fontsize=10)
plt.ylabel(r'$R_{in situ}$ []', fontsize=10)
plt.savefig('../figures/ratioR_latitude.pdf', format='pdf')
plt.show()