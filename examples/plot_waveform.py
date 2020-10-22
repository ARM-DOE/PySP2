import pysp2
import time
import matplotlib.pyplot as plt
import numpy as np
import os

if __name__ == "__main__":
    quicklook_output_dir = 'C:\\Users\\rjackson\\Documents\\SP2Data4ANL\\AquaDag_PSL_cal\\quicklooks\\'
    ptime = time.time()
    my_binary = pysp2.io.read_sp2('C:\\Users\\rjackson\\Documents\\SP2Data4ANL\\AquaDag_PSL_cal\\20190703\\20190703x006.sp2b',
                                  debug=True)
    print("## SP2 record read in " + str(time.time()-ptime) + "s")

    my_config = pysp2.io.read_config('C:\\Users\\rjackson\Documents\\SP2Data4ANL\\AquaDag_PSL_cal\\20190703\\20190703172341.ini')
    print(my_config['Acquisition']['Pre-Trig Points'])

    ptime = time.time()
    my_binary = pysp2.util.gaussian_fit(my_binary, my_config, parallel=False)
    for i in range(len(my_binary['Base_ch0'].values)):
        my_record = i
        display = pysp2.vis.plot_wave(my_binary, my_record, 0)
        base = my_binary['Base_ch0'].values[i]
        amp = my_binary['FtAmp_ch0'].values[i]
        display.axes[0].set_title("Base_ch0 = %5.2f FtAmp_ch0 = %5.2f" % (base, amp))
        display.fig.savefig(os.path.join(quicklook_output_dir, '%d.png' % i))
        plt.close(display.fig)
