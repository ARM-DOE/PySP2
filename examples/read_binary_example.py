import pysp2
import time
import matplotlib.pyplot as plt
import numpy as np

from glob import glob
from distributed import LocalCluster, Client

if __name__ == "__main__":
    quicklook_output_dir = '/nfs/gce/projects/digr/SP2Data4ANL/Pysp2_test/'
    ptime = time.time()
    for i in range(1,16):
        my_binary = pysp2.io.read_sp2('/nfs/gce/projects/digr/SP2Data4ANL/20181110/20181110x%03d.sp2b' % i,
                                      debug=True)
        print("## SP2 record read in " + str(time.time()-ptime) + "s")

        my_config = pysp2.io.read_config('/nfs/gce/projects/digr/SP2Data4ANL/20181110/20181110114046.ini')
        print(my_config['Acquisition']['Pre-Trig Points'])

        ptime = time.time()
        my_binary = pysp2.util.gaussian_fit(my_binary, my_config, parallel=True)
        print(my_binary)
        print("## SP2 fits made in %4.2f s" % (time.time()-ptime))
        pysp2.io.write_dat(my_binary, '/nfs/gce/projects/digr/SP2Data4ANL/pyprocessing/20181110x%03d.dat' % i )

    # for i in range(100):
    #      coeffs = [my_binary['FtAmp_ch0'].values[i],
    #                my_binary['FtPos_ch0'].values[i],
    #                my_binary['PkFWHM_ch0'].values[i]/(2.35482/np.sqrt(2)),
    #                my_binary['Base_ch0'].values[i]]
    #      plt.plot(my_binary['Data_ch0'].values[i, :], label='Obs', color='k')
    #      plt.plot(pysp2.util._gaus(np.arange(0, 100), *coeffs), label='fit', color='b')
    #      plt.legend()
    #      plt.title('Amplitude = %3.2f Peak = %3.2f Width=%3.2f Base=%3.2f' %
    #                (coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
    #      plt.savefig(quicklook_output_dir + '/ch0/' + str(i) + '.png')
    #      plt.close()
    #
    #      coeffs = [my_binary['FtAmp_ch4'].values[i],
    #                my_binary['FtPos_ch4'].values[i],
    #                my_binary['PkFWHM_ch4'].values[i]/(2.35482/np.sqrt(2)),
    #                my_binary['Base_ch4'].values[i]]
    #      plt.plot(my_binary['Data_ch4'].values[i, :], label='Obs', color='k')
    #      plt.plot(pysp2.util._gaus(np.arange(0, 100), *coeffs), label='fit', color='b')
    #      plt.legend()
    #      plt.title('Amplitude = %3.2f Peak = %3.2f Width=%3.2f Base=%3.2f' %
    #                (coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
    #      plt.savefig(quicklook_output_dir + '/ch4/' + str(i) + '.png')
    #      plt.close()
    #
    #      coeffs = [my_binary['PkHt_ch1'].values[i],
    #                my_binary['PkSplitPos_ch1'].values[i],
    #                my_binary['PkStart_ch1'].values[i],
    #                my_binary['PkEnd_ch1'].values[i],
    #                my_binary['Base_ch1'].values[i]]
    #      plt.plot(my_binary['Data_ch1'].values[i, :], label='Obs', color='k')
    #      plt.legend()
    #      plt.title('Height = %3.2f Peak = %3.2f Start=%3.2f End=%3.2f Base=%3.2f' %
    #                (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]))
    #      plt.savefig(quicklook_output_dir + '/ch1/' + str(i) + '.png')
    #      plt.close()
    #      coeffs = [my_binary['PkHt_ch2'].values[i],
    #                my_binary['PkSplitPos_ch2'].values[i],
    #                my_binary['PkStart_ch2'].values[i],
    #                my_binary['PkEnd_ch2'].values[i],
    #                my_binary['Base_ch2'].values[i]]
    #      plt.plot(my_binary['Data_ch2'].values[i, :], label='Obs', color='k')
    #      plt.legend()
    #      plt.title('Height = %3.2f Peak = %3.2f Start=%3.2f End=%3.2f Base=%3.2f' %
    #                (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]))
    #      plt.savefig(quicklook_output_dir + '/ch2/' + str(i) + '.png')
    #      plt.close()
    #
    #      coeffs = [my_binary['PkHt_ch3'].values[i],
    #                my_binary['PkPos_ch3'].values[i],
    #                my_binary['PkSplitPos_ch3'].values[i],
    #                my_binary['Base_ch3'].values[i]]
    #      plt.plot(my_binary['Data_ch3'].values[i, :], label='Obs', color='k')
    #      plt.legend()
    #      plt.title('Height = %3.2f Peak = %3.2f Split=%3.2f Base=%3.2f' %
    #                (coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
    #      plt.savefig(quicklook_output_dir + '/ch3/' + str(i) + '.png')
    #      plt.close()
    #
    #      coeffs = [my_binary['PkHt_ch5'].values[i],
    #                my_binary['PkSplitPos_ch5'].values[i],
    #                my_binary['PkStart_ch5'].values[i],
    #                my_binary['PkEnd_ch5'].values[i],
    #                my_binary['Base_ch5'].values[i]]
    #      plt.plot(my_binary['Data_ch5'].values[i, :], label='Obs', color='k')
    #      plt.legend()
    #      plt.title('Height = %3.2f Peak = %3.2f Start=%3.2f End=%3.2f Base=%3.2f' %
    #                (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]))
    #      plt.savefig(quicklook_output_dir + '/ch5/' + str(i) + '.png')
    #      plt.close()
    #      coeffs = [my_binary['PkHt_ch6'].values[i],
    #                my_binary['PkSplitPos_ch6'].values[i],
    #                my_binary['PkStart_ch6'].values[i],
    #                my_binary['PkEnd_ch6'].values[i],
    #                my_binary['Base_ch6'].values[i]]
    #      plt.plot(my_binary['Data_ch6'].values[i, :], label='Obs', color='k')
    #      plt.legend()
    #
    #      plt.title('Height = %3.2f Peak = %3.2f Start=%3.2f End=%3.2f Base=%3.2f' %
    #                (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]))
    #      plt.savefig(quicklook_output_dir + '/ch6/' + str(i) + '.png')
    #      plt.close()
    #
    #      coeffs = [my_binary['PkHt_ch7'].values[i],
    #                my_binary['PkPos_ch7'].values[i],
    #                my_binary['PkSplitPos_ch7'].values[i],
    #                my_binary['Base_ch7'].values[i]]
    #      plt.plot(my_binary['Data_ch7'].values[i, :], label='Obs', color='k')
    #      plt.legend()
    #      plt.title('Height = %3.2f Peak = %3.2f Split=%3.2f Base=%3.2f' %
    #                (coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
    #      plt.savefig(quicklook_output_dir + '/ch7/' + str(i) + '.png')
    #      plt.close()
    #
