import pysp2
import pandas as pd
import glob
import numpy as np

dat_path = '/lustre/or-hydra/cades-arm/rjackson/nsaaossp2X2.b1/psd/20200218/*.dat'
header = open('header_arm_dat_nsa', 'r')
#out_file = open('sp2-nsa-02-18-2020-60s.dat', 'w')
#print(glob.glob(dat_path, recursive=True))
ds = pysp2.io.read_arm_dat(dat_path)

while True:
    line = header.readline()
    if not line:
        break
    out_file.write(line)

def format_time(x):
    return '%11.1f' % x

#ds['SP2_Dmin'][200:] = np.nan
#ds['SP2_Dgeo'][200:] = np.nan
#ds['SP2_Dmax'][200:] = np.nan
ds = ds.fillna('')
print(ds)
#ds["SP2_datetime_in_sec"] = ds["SP2_datetime_in_sec"].apply(format_time)
#ds = ds.sort_values("SP2_datetime_in_sec")
#ds = ds.to_dataframe()
ds.to_csv(out_file, float_format = "%.8g", encoding='utf-8', sep='\t')
out_file.close()
