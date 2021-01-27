import pysp2
import pandas as pd

dat_path = '/home/rjackson/1s_data_sample/*.dat'
ds = pysp2.io.read_dat(dat_path, 'conc')
print(ds)
ds = ds.to_dataframe()
ds.to_csv('sp2-nsa-2020-1s.dat', index=False, float_format = "%.8g", encoding='utf-8', sep='\t')
