import pysp2
import time
import numpy as np

if __name__ == "__main__":
    ptime = time.time()
    for i in range(1, 16):
        my_binary = pysp2.io.read_sp2('/nfs/gce/projects/digr/SP2Data4ANL/20181110/20181110x%03d.sp2b' % i,
                                      debug=True)
        print("## SP2 record read in " + str(time.time()-ptime) + "s")
        my_config = pysp2.io.read_config('/nfs/gce/projects/digr/SP2Data4ANL/20181110/20181110114046.ini')
        ptime = time.time()
        my_binary = pysp2.util.gaussian_fit(my_binary, my_config, parallel=False)
        print("## SP2 fits made in %4.2f s" % (time.time()-ptime))
        pysp2.io.write_dat(my_binary, '/nfs/gce/projects/digr/SP2Data4ANL/pyprocessing/20181110x%03d.dat' % i )

