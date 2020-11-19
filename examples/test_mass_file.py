import xarray as xr
import matplotlib.pyplot as plt
import pysp2
masses = xr.open_dataset('SP2_masses.nc')
print(masses)
pysp2.io.write_dat_concs(masses, 'SP2test.dat')
#masses['NumConcIncanScat'].plot()
#plt.show()
masses.close()