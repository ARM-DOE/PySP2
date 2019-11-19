import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np

my_dat_files = sorted(glob.glob('/nfs/gce/projects/digr/SP2Data4ANL/pyprocessing/20181110x0*.dat'))
igor_dat_files = sorted(glob.glob('/nfs/gce/projects/digr/SP2Data4ANL/20181110/20181110x*.dat'))
igor_dat_files = sorted(list(set(igor_dat_files) - set(my_dat_files)))

python_ds = pd.concat([pd.read_csv(x, sep='\t', skiprows=2) for x in my_dat_files])
igor_ds = pd.concat([pd.read_csv(x, sep='\t', skiprows=2) for x in igor_dat_files])

print(python_ds.columns.values)
variable = 'FtAmp_ch0'
plt.figure(figsize=(10,10))
plt.scatter(python_ds[variable], igor_ds[variable])
minv = python_ds[variable].min()
maxv = python_ds[variable].max()
plt.plot([-1e2, 1e9], [-1e9, 1e9], linewidth=2, color='k')
plt.xlim([minv, maxv])
plt.ylim([minv, maxv])
plt.xlabel(variable + " Python code")
plt.ylabel(variable + " Igor code")
plt.show()

diffs = abs(python_ds[variable].values-igor_ds[variable].values)
how_many_diff = len(np.where(diffs > 10)[0])

plt.figure(figsize=(10,10))
plt.plot(python_ds['TimeWave'].values, python_ds['FtAmp_ch0'].values)
plt.plot(python_ds['TimeWave'].values, igor_ds['FtAmp_ch0'].values)
plt.ylim([0, 60000])
plt.show()