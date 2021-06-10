import pysp2
import pandas as pd
import glob
import numpy as np
import sys

year = sys.argv[1]
day = sys.argv[2]
site = sys.argv[3]

if site == 'nsa':
    datastream = 'nsaaossp2X2.b1'
elif site == 'mos':
    datastream = 'mosaossp2M1.b1'

dat_path = '/lustre/or-scratch/cades-arm/rjackson/%s/final_psds/*%s%s*.dat' % (datastream, year, day)
header = open('header_arm_dat_%s' % site, 'r')
out_file = open('sp2-%s-%s-%s-60s.dat' % (site, day, year), 'w')
file_list = sorted(glob.glob(dat_path, recursive=True))
ds_list = []
def format_time(x):
    return '%11.1f' % x

for fname in file_list:
    ds_temp = pysp2.io.read_arm_dat(fname)
    #ds_temp.index = ds_temp.index.map(format_time)
    ds_list.append(ds_temp)
    

ds = pd.concat(ds_list)
print(len(ds_list))

while True:
    line = header.readline()
    if not line:
        break
    out_file.write(line)

SP2_min = np.arange(7.5, 1007.5, 5.)
SP2_geo = np.arange(10., 1010., 5.)
SP2_max = np.arange(12.5, 1012.5, 5.)
ds["SP2_Dmin"].values[0:len(SP2_min)] = SP2_min
ds["SP2_Dgeo"].values[0:len(SP2_geo)] = SP2_geo
ds["SP2_Dmax"].values[0:len(SP2_max)] = SP2_max
ds['SP2_Dmin'].values[200:] = np.nan
ds['SP2_Dgeo'].values[200:] = np.nan
ds['SP2_Dmax'].values[200:] = np.nan
ds = ds.fillna('')
ds = ds.sort_index()
print(ds.index[0:10])
#ds = ds.to_dataframe()
ds.to_csv(out_file, float_format = "%.10g", encoding='utf-8', sep='\t')
out_file.close()
