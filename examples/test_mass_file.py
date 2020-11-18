import xarray as xr
import matplotlib.pyplot as plt

masses = xr.open_dataset('SP2_masses.nc')
print(masses)
masses['ScatMassEnsemble'].plot()
plt.show()
masses.close()