import xarray as xr
import matplotlib.pyplot as plt
import pysp2
import numpy as np

dat_py = pysp2.io.read_dat('/nfs/gce/projects/digr/SP2Data4ANL/pyprocessing/*.dat', type='particle')
#dat_py.load()
dat_py = pysp2.util.calc_diams_masses(dat_py)
dat_py.to_netcdf('PyProcessing_properties.nc')

dat_igor = pysp2.io.read_dat('/nfs/gce/projects/digr/SP2Data4ANL/20181110/20181110x*.dat', type='particle')
#dat_igor.load()
dat_igor = pysp2.util.calc_diams_masses(dat_igor)
dat_igor.to_netcdf('IgorProcessing_properties.nc')

hist, bins = np.histogram(dat_py['sootDiam'], bins=np.linspace(0, 500, 100), normed=True)
hist_igor, bins = np.histogram(dat_igor['sootDiam'], bins=np.linspace(0, 500, 100), normed=True)

fig, ax = plt.subplots(1, 1)
ax.step(bins[1:], hist*bins[1:], label='Python')
ax.step(bins[1:], hist_igor*bins[1:], label='Igor')
ax.set_xlabel('Incadesence diameter [nm]')
ax.set_ylabel('p.d.f [nm]')
ax.legend()
plt.show(fig)
hist, bins = np.histogram(dat_py['ScatDiaSO4'], bins=np.linspace(0, 700, 100), normed=True)
hist_igor, bins = np.histogram(dat_igor['ScatDiaSO4'], bins=np.linspace(0, 700, 100), normed=True)

fig, ax = plt.subplots(1, 1)
ax.step(bins[1:], hist, label='Python')
ax.step(bins[1:], hist_igor, label='Igor')
ax.set_xlabel('scattering diameter [nm]')
ax.set_ylabel('p.d.f')
ax.legend()
plt.show(fig)


