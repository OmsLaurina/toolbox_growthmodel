"""

Plot the in situ biomasses. Create the figure 2 of the paper and the left panel of the figure 9.

"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker
import scienceplots

plt.style.use(['science','no-latex'])
plt.close('all')

# Load data from outputs_insitu.py code
pico = np.loadtxt('../outputs/pico.txt')
nano = np.loadtxt('../outputs/nano.txt')
micro = np.loadtxt('../outputs/micro.txt')
lat = np.loadtxt('../outputs/lat.txt')
lon = np.loadtxt('../outputs/lon.txt')
f_BioM = np.loadtxt('../outputs/f_biomass.txt')
f_biomass_syne,f_biomass_nano,f_biomass_micro = f_BioM.T

# Figure 1 (figure 2 of the paper): Biomasses distribution
fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(1, 3, 1) # PICO
ax2 = fig.add_subplot(1, 3, 2) # MICRO
ax3 = fig.add_subplot(1, 3, 3) # NANO

scatter1 = ax1.scatter(lon, lat, c=pico, cmap='viridis')
ax1.set_title(r'$PICO$' + '\n$f_{BioM}=' + str(round(f_biomass_syne, 1)) + '$')
ax1.set_xlabel('°E Longitude')
ax1.set_ylabel('°N Latitude')
ax1.set_aspect('equal')
cbar1 = plt.colorbar(scatter1, ax=ax1, label='Masses [mmolC.m$^{-3}$]')
ax1.axhline(y=38.5, color='lightcoral', linestyle='-',linewidth=1.5)
tick_locator = ticker.MaxNLocator(nbins=10)
cbar1.locator = tick_locator
cbar1.update_ticks()

scatter2 = ax2.scatter(lon, lat, c=micro, cmap='viridis')
ax2.set_title(r'$MICRO$' + '\n$f_{BioM}=' + str(round(f_biomass_micro, 1)) + '$')
ax2.set_xlabel('°E Longitude')
ax2.set_ylabel('°N Latitude')
ax2.set_aspect('equal')
cbar2 = plt.colorbar(scatter2, ax=ax2, label='Masses [mmolC.m$^{-3}$]')
ax2.axhline(y=38.5, color='lightcoral', linestyle='-',linewidth=1.5)
tick_locator = ticker.MaxNLocator(nbins=10)
cbar2.locator = tick_locator
cbar2.update_ticks()

scatter3 = ax3.scatter(lon, lat, c=nano, cmap='viridis')
ax3.set_title(r'$NANO$' + '\n$f_{BioM}=' + str(round(f_biomass_nano, 1)) + '$')
ax3.set_xlabel('°E Longitude')
ax3.set_ylabel('°N Latitude')
ax3.set_aspect('equal')
cbar3 = plt.colorbar(scatter3, ax=ax3, label='Masses [mmolC.m$^{-3}$]')
ax3.axhline(y=38.5, color='lightcoral', linestyle='-',linewidth=1.5)
tick_locator = ticker.MaxNLocator(nbins=10)
cbar3.locator = tick_locator
cbar3.update_ticks()

plt.tight_layout()
plt.show()
plt.savefig('../figures/insitu_biomasses.pdf', format='pdf')

# Calculate the polynomial regression of the in situ R ratio (pico/mean(to))
tot = np.loadtxt('../outputs/tot.txt')
degree = 3
coefficients = np.polyfit(lat, pico/np.mean(tot), degree)
polynomial = np.poly1d(coefficients)
x_fit = np.linspace(min(lat), max(lat), 1000)
y_fit = polynomial(x_fit)

# Figure 2 (figure 9 of the paper): R ratio as a function of the latitude
plt.figure(2)
plt.scatter(lat, pico/np.mean(tot), c=tot, cmap='Wistia')
plt.plot(x_fit, y_fit, 'k-', label=f'Polynomial Fit (Degree {degree})')
plt.axvline(x=38.5, color='red', linestyle='--',linewidth=2)
plt.grid()
plt.colorbar(label=r'${TotBioM}$ [mmolCm$^{-3}$]')
plt.xlabel('°N Latitude', fontsize=10)
plt.ylabel(r'$R_{in situ}$ []', fontsize=10)
plt.savefig('../figures/ratioR_latitude.pdf', format='pdf')