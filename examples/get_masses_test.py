import pysp2
import matplotlib.pyplot as plt

hk_file = 'nsaaossp2auxX2.00.20200220.001001.raw.20200219000000.hk'
dat_file = 'nsaaossp2X2.00.20200219.150801.raw.20200219x018.sp2b.dat'
config_file = 'nsaaossp2auxX2.00.20200219.154225.raw.20200116000000.ini'

hk_ds = pysp2.io.read_hk_file(hk_file)
ds = pysp2.io.read_dat(dat_file, type='particle')
config = pysp2.io.read_config(config_file)
particles = pysp2.util.calc_diams_masses(ds)
out_mass = pysp2.util.process_psds(particles, hk_ds, config)
out_mass.to_netcdf('SP2_masses.nc')
