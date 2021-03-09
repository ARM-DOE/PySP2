import pysp2
import os
import matplotlib.pyplot as plt
import dask.bag as db

from glob import glob
from dask_jobqueue import SLURMCluster
from distributed import Client, wait
from datetime import datetime, timedelta

sp2b_path = '/lcrc/group/earthscience/rjackson/sp2_sample/'
sp2baux_path = '/lcrc/group/earthscience/rjackson/sp2_sample/'

dat_path = '/lcrc/group/earthscience/rjackson/sp2_sample/particle_info/'
out_path = '/lcrc/group/earthscience/rjackson/sp2_sample/psds/'
quicklook_path = '/lcrc/group/earthscience/rjackson/psds/quicklooks/'

def process_file(fname, hk_file, ini_file, my_day):
    base = os.path.split(fname)[-1]
    #if os.path.isfile(out_path + '/' + my_day + '/' + base + '.1s.dat'):
    #    continue
    if not os.path.isdir(out_path + '/' + my_day):
        os.makedirs(out_path + '/' + my_day)

    if not os.path.isdir(quicklook_path + '/' + my_day):
        os.makedirs(quicklook_path + '/' + my_day)
    print("Processing %s" % fname)
    my_hk = pysp2.io.read_hk_file(hk_file)
    my_config = pysp2.io.read_config(ini_file)
    my_binary = pysp2.io.read_dat(fname, type='particle')
    particles = pysp2.util.calc_diams_masses(my_binary, factor=1)
    out_mass = pysp2.util.process_psds(
        particles, my_hk, my_config, avg_interval=1)
    pysp2.io.write_dat_concs(out_mass,
        out_path + '/' + my_day + '/' + base + '.1s.dat')

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
    my_binary.close()
    del particles, my_binary, out_mass


def process_day(my_day):
    all_files = sorted(glob(dat_path + '*' + my_day + '*.dat'))
    ini_file = sorted(glob(sp2baux_path + '*' + my_day + '*.ini'))[0]
    
    my_day_plus_one = datetime.strptime(my_day, '%Y%m%d') + timedelta(days=1)
    my_day_plus_one = my_day_plus_one.strftime('%Y%m%d')
    hk_files = sorted(glob(sp2baux_path + '*.raw.' + my_day + '*.hk'))
    my_func = lambda x: process_file(x, hk_files[0], ini_file, my_day) 
    file_bag = db.from_sequence(all_files)
    file_bag.map(my_func).compute()

# Get all of the unique dates
all_sp2_files = glob(sp2b_path + '*.sp2b')
#sp2_date_list = [x.split(".")[3] for x in all_sp2_files]
#sp2_date_list = sorted(list(set(sp2_date_list)))
sp2_date_list = ['20200218']
#process_day('20200218')
#print(sp2_date_list)
cluster = SLURMCluster(processes=9, cores=32, walltime='1:00:00', memory='128GB', name='dask-worker', project='rainfall')
cluster.scale(18)
client = Client(cluster)

print("Waiting for workers before starting processing...")
client.wait_for_workers(9)
print(client)
process_day('20200218')
#results = client.map(process_day, sp2_date_list)
#wait(results)
del client
