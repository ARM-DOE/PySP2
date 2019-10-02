import pysp2
import time
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    quicklook_output_dir = '/nfs/gce/projects/digr/SP2Data4ANL/Pysp2_test/'
    ptime = time.time()
    my_binary = pysp2.io.read_sp2('/nfs/gce/projects/digr/SP2Data4ANL/20181110/20181110x001.sp2b',
                                  debug=True)
    print("## SP2 record read in " + str(time.time()-ptime) + "s")

    my_config = pysp2.io.read_config('/nfs/gce/projects/digr/SP2Data4ANL/20181110/20181110114046.ini')
    print(my_config['Acquisition']['Pre-Trig Points'])

    ptime = time.time()
    #my_binary = pysp2.util.gaussian_fit(my_binary, my_config, parallel=True)
    my_record = 3500
    pysp2.vis.plot_wave(my_binary, my_record, 1)
    plt.show()
