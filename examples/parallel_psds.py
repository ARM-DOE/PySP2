import pysp2
import xarray as xr
import os
import matplotlib.pyplot as plt
import dask.bag as db
import sys

from glob import glob
from dask_jobqueue import PBSCluster
from distributed import Client, wait
from datetime import datetime, timedelta

sp2b_path = '/lustre/or-hydra/cades-arm/rjackson/nsaaossp2X2.00/'
sp2baux_path = '/lustre/or-hydra/cades-arm/rjackson/nsaaossp2auxX2.00/'

dat_path = '/lustre/or-hydra/cades-arm/rjackson/nsaaossp2X2.a1/'
out_path = '/lustre/or-hydra/cades-arm/rjackson/nsaaossp2X2.b1/psd/'
quicklook_path = '/lustre/or-hydra/cades-arm/rjackson/nsaaossp2X2.b1/quicklooks/'
cal_path = '/home/rjackson/PySP2/examples/Unit25CAL_20201016_PostDMT_AQ_FullereneEquiv.txt'


def process_file(fname, hk_file, ini_file, my_day):
    base = os.path.split(fname)[-1]
    #if os.path.isfile(out_path + '/' + my_day + '/' + base + '.1s.dat'):
    #    continue
    if not os.path.isdir(out_path + '/' + my_day):
        os.makedirs(out_path + '/' + my_day)

    if not os.path.isdir(quicklook_path + '/' + my_day):
        os.makedirs(quicklook_path + '/' + my_day)
    print("Processing %s" % fname)
    my_binary = pysp2.io.read_dat(fname, type='particle')
    Globals = pysp2.util.DMTGlobals(cal_path)
    particles = pysp2.util.calc_diams_masses(
        my_binary, factor=1.3, Globals=Globals)
    return particles
    
def process_day(my_day):
    all_files = sorted(glob(dat_path + '*.raw.' + my_day + 'x*.dat'))
    ini_file = sorted(glob(sp2baux_path + '*' + my_day + '*.ini'))[0]
    
    my_day_plus_one = datetime.strptime(my_day, '%Y%m%d') + timedelta(days=1)
    my_day_plus_one = my_day_plus_one.strftime('%Y%m%d')
    hk_files = sorted(glob(sp2baux_path + '*raw.' + my_day + '*.hk'))
    if len(hk_files) == 0:
        return
    try:
        my_func = lambda x: process_file(x, hk_files, ini_file, my_day) 
    except IndexError:
        return
    for i in range(0, len(all_files), 10):
        next_i = min([i + 10, len(all_files)])
        file_list = all_files[i:next_i]
        base = os.path.split(file_list[0])[-1] 
        particles = []
        for f in file_list:
            particles.append(my_func(f))
        particles = xr.concat(particles, dim='index')
        
        my_hk = []
        for f in hk_files:
            my_hk.append(pysp2.io.read_hk_file(f))
        my_hk = xr.concat(my_hk, dim='time') 
        my_config = pysp2.io.read_config(ini_file)

        out_mass = pysp2.util.process_psds(
            particles, my_hk, my_config, avg_interval=60)

        fig, ax = plt.subplots(2, 1, figsize=(10, 7))
        out_mass.NumConcScat.plot(ax=ax[0], linewidth=2, color='r',
                                  label='scattering')
        out_mass.NumConcIncan.plot(ax=ax[0], linewidth=2, color='b',
                                  label='incandescence')
        ax[0].legend()
        ax[0].set_ylim([0, 50])
        ax[0].set_ylabel('Concentration [cm-3]')

        out_mass.MassScat2.plot(ax=ax[1], linewidth=2, color='r',
                                label='scattering')
        out_mass.MassIncand2total.plot(ax=ax[1], linewidth=2, color='b',
                                       label='incandescence')
        ax[1].set_ylim([0, 300])
        ax[1].legend()
        ax[1].set_ylabel('Mass concentration [ng m-3]')
        fig.savefig(
            quicklook_path + '/' + my_day + '/' + base + 'concsmasses.png')
        pysp2.io.write_dat_concs_arm(out_mass,
            out_path + '/' + my_day + '/' + base + '.60s.dat',
            location="Point Barrow", lat_lon_string="71.2906 N 156.7886 W")
        del particles
        plt.close(fig)
    


if __name__ == '__main__':
    # Get all of the unique dates
    all_sp2_files = glob(sp2b_path + '*' + sys.argv[1] + '*.sp2b')
    sp2_date_list = [x.split(".")[3] for x in all_sp2_files]
    sp2_date_list = sorted(list(set(sp2_date_list)))
    
    cluster = PBSCluster(processes=9, cores=36, walltime='4:00:00',
                         memory='270GB', name='dask-worker',
                         queue='arm_high_mem',
                         project='arm', job_extra=['-W group_list=cades-arm'],
                         interface='ib0',
                         extra=['--no-dashboard'])
    cluster.scale(72)

    client = Client(cluster)

    print("Waiting for workers before starting processing...")
    client.wait_for_workers(9)
    #for day in sp2_date_list:
    #    print("Processing %s" % day)
    #    process_day(day)
    #print(sp2_date_list)
    #process_day('20200218')
    results = client.map(process_day, sp2_date_list)
    wait(results)
    del client
