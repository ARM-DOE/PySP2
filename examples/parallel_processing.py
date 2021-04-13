import pysp2
import os

from glob import glob
from dask_jobqueue import PBSCluster
from distributed import Client, wait

sp2b_path = '/lustre/or-hydra/cades-arm/rjackson/nsaaossp2X2.00/'
sp2baux_path = '/lustre/or-hydra/cades-arm/rjackson/nsaaossp2auxX2.00/'

out_path = '/lustre/or-hydra/cades-arm/rjackson/nsaaossp2X2.a1/'

def process_day(my_day):
    all_files = sorted(glob(sp2b_path + '*' + my_day + '*.sp2b'))
    ini_file = sorted(glob(sp2baux_path + '*' + my_day + '*.ini'))[0]
    
    for fname in all_files:
        base = os.path.split(fname)[-1] 
        if os.path.isfile(out_path + base + '.particle.dat'):
            continue
        print("Processing %s" % fname)
        my_binary = pysp2.io.read_sp2(fname)
        my_config = pysp2.io.read_config(ini_file)
        my_binary = pysp2.util.gaussian_fit(
            my_binary, my_config, parallel=False)
        pysp2.io.write_dat(my_binary, out_path + base + '.particle.dat')

# Get all of the unique dates
all_sp2_files = glob(sp2b_path + '*.sp2b')
sp2_date_list = [x.split(".")[3] for x in all_sp2_files]
sp2_date_list = sorted(list(set(sp2_date_list)))
#sp2_date_list = ['20200218']
#process_day('20200218')
print(sp2_date_list)
cluster = PBSCluster(processes=6, cores=36, walltime='5:00:00', 
                     memory='270GB', name='dask-worker', queue='arm_high_mem',
                     project='arm', job_extra=['-W group_list=cades-arm'],
                     interface='ib0', 
                     extra=['--no-dashboard'])
cluster.scale(36*6)
client = Client(cluster)

print("Waiting for workers before starting processing...")
client.wait_for_workers(9)
print(client)
results = client.map(process_day, sp2_date_list)
wait(results)
#del client
 

     

