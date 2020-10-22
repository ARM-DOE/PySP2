import pysp2
import time
import numpy as np

if __name__ == "__main__":
    ptime = time.time()
    for i in range(1, 21):
        my_binary = pysp2.io.read_sp2('/nfs/gce/projects/digr/SP2Data4ANL/AquaDag_PSL_cal/20190703/20190703x%03d.sp2b' % i,
                                      debug=True)
        print("## SP2 record read in " + str(time.time()-ptime) + "s")
        my_config = pysp2.io.read_config('/nfs/gce/projects/digr/SP2Data4ANL/AquaDag_PSL_cal/20190703/20190703172341.ini')
        ptime = time.time()
        my_binary = pysp2.util.gaussian_fit(my_binary, my_config, parallel=False)
        print("## SP2 fits made in %4.2f s" % (time.time()-ptime))
        pysp2.io.write_dat(my_binary, '/nfs/gce/projects/digr/SP2Data4ANL/AquaDag_PSL_cal/Py_generated_numbers/20190703x%03d.dat' % i)

