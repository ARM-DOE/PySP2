import pysp2
import pandas as pd

dat_path = '/lcrc/group/earthscience/rjackson/sp2_sample/psds/20200218/*.dat'
ds = pysp2.io.read_dat(dat_path, 'conc')
print(ds)
ds = ds.to_dataframe()
ds.to_csv('sp2-nsa-2020-60s.dat', index=False, float_format = "%.8g", encoding='utf-8', sep='\t')
