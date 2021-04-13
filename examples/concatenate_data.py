import pysp2
import pandas as pd
import glob
import numpy as np

dat_path = '/lustre/or-hydra/cades-arm/rjackson/mosaossp2M1.b1/psd/20191113/*.dat'
header = open('header_arm_dat_mos', 'r')
out_file = open('sp2-mosaic-11-13-2020-60s.dat', 'w')
#print(glob.glob(dat_path, recursive=True))
ds = pysp2.io.read_arm_dat(dat_path)

while True:
    line = header.readline()
    if not line:
        break
    out_file.write(line)

#ds['SP2_Dmin'][200:] = np.nan
#ds['SP2_Dgeo'][200:] = np.nan
#ds['SP2_Dmax'][200:] = np.nan
ds = ds.fillna('')
print(ds)
#ds = ds.to_dataframe()
ds.to_csv(out_file, index=False, float_format = "%.8g", encoding='utf-8', sep='\t')
out_file.close()
