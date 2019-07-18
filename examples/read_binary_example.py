import pysp2
import matplotlib.pyplot as plt
import numpy as np

quicklook_output_dir = '/home/rjackson/Downloads/SP2Data4ANL/Pysp2_test/'
my_binary = pysp2.io.read_sp2('/home/rjackson/Downloads/SP2Data4ANL/20190518x001.sp2b',
                              debug=True)
my_config = pysp2.io.read_config('/home/rjackson/Downloads/SP2Data4ANL/20190518000000.ini')
print(my_config['Acquisition']['Pre-Trig Points'])


my_binary = pysp2.util.gaussian_fit(my_binary, my_config, num_records=300)


for i in range(300):
    coeffs = [my_binary['amplitude_ch0'].values[i],
              my_binary['peak_pos_ch0'].values[i],
              my_binary['width_ch0'].values[i],
              my_binary['base_ch0'].values[i]]
    plt.plot(my_binary['Data_ch0'].values[i, :], label='Obs', color='k')
    plt.plot(pysp2.util._gaus(np.arange(0, 100), *coeffs), label='fit', color='b')
    plt.legend()
    plt.title('Amplitude = %3.2f Peak = %3.2f Width=%3.2f Base=%3.2f' %
              (coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
    plt.savefig(quicklook_output_dir + '/ch0/' + str(i) + '.png', dpi=300)
    plt.close()

    coeffs = [my_binary['amplitude_ch4'].values[i],
              my_binary['peak_pos_ch4'].values[i],
              my_binary['width_ch4'].values[i],
              my_binary['base_ch4'].values[i]]
    plt.plot(my_binary['Data_ch4'].values[i, :], label='Obs', color='k')
    plt.plot(pysp2.util._gaus(np.arange(0, 100), *coeffs), label='fit', color='b')
    plt.legend()
    plt.title('Amplitude = %3.2f Peak = %3.2f Width=%3.2f Base=%3.2f' %
              (coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
    plt.savefig(quicklook_output_dir + '/ch4/' + str(i) + '.png', dpi=300)
    plt.close()

    coeffs = [my_binary['height_ch1'].values[i],
              my_binary['peak_pos_ch1'].values[i],
              my_binary['peak_start_ch1'].values[i],
              my_binary['peak_end_ch1'].values[i],
              my_binary['base_ch1'].values[i]]
    plt.plot(my_binary['Data_ch1'].values[i, :], label='Obs', color='k')
    plt.legend()
    plt.title('Height = %3.2f Peak = %3.2f Start=%3.2f End=%3.2f Base=%3.2f' %
              (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]))
    plt.savefig(quicklook_output_dir + '/ch1/' + str(i) + '.png', dpi=300)
    plt.close()
    coeffs = [my_binary['height_ch2'].values[i],
              my_binary['peak_pos_ch2'].values[i],
              my_binary['peak_start_ch2'].values[i],
              my_binary['peak_end_ch2'].values[i],
              my_binary['base_ch2'].values[i]]
    plt.plot(my_binary['Data_ch2'].values[i, :], label='Obs', color='k')
    plt.legend()
    plt.title('Height = %3.2f Peak = %3.2f Start=%3.2f End=%3.2f Base=%3.2f' %
              (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]))
    plt.savefig(quicklook_output_dir + '/ch2/' + str(i) + '.png', dpi=300)
    plt.close()

    coeffs = [my_binary['height_ch3'].values[i],
              my_binary['peak_pos_ch3'].values[i],
              my_binary['peak_split_ch3'].values[i],
              my_binary['base_ch3'].values[i]]
    plt.plot(my_binary['Data_ch3'].values[i, :], label='Obs', color='k')
    plt.legend()
    plt.title('Height = %3.2f Peak = %3.2f Split=%3.2f Base=%3.2f' %
              (coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
    plt.savefig(quicklook_output_dir + '/ch3/' + str(i) + '.png', dpi=300)
    plt.close()

    coeffs = [my_binary['height_ch5'].values[i],
              my_binary['peak_pos_ch5'].values[i],
              my_binary['peak_start_ch5'].values[i],
              my_binary['peak_end_ch5'].values[i],
              my_binary['base_ch5'].values[i]]
    plt.plot(my_binary['Data_ch5'].values[i, :], label='Obs', color='k')
    plt.legend()
    plt.title('Height = %3.2f Peak = %3.2f Start=%3.2f End=%3.2f Base=%3.2f' %
              (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]))
    plt.savefig(quicklook_output_dir + '/ch5/' + str(i) + '.png', dpi=300)
    plt.close()
    coeffs = [my_binary['height_ch6'].values[i],
              my_binary['peak_pos_ch6'].values[i],
              my_binary['peak_start_ch6'].values[i],
              my_binary['peak_end_ch6'].values[i],
              my_binary['base_ch6'].values[i]]
    plt.plot(my_binary['Data_ch6'].values[i, :], label='Obs', color='k')
    plt.legend()

    plt.title('Height = %3.2f Peak = %3.2f Start=%3.2f End=%3.2f Base=%3.2f' %
              (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]))
    plt.savefig(quicklook_output_dir + '/ch6/' + str(i) + '.png', dpi=300)
    plt.close()

    coeffs = [my_binary['height_ch7'].values[i],
              my_binary['peak_pos_ch7'].values[i],
              my_binary['peak_split_ch7'].values[i],
              my_binary['base_ch7'].values[i]]
    plt.plot(my_binary['Data_ch7'].values[i, :], label='Obs', color='k')
    plt.legend()
    plt.title('Height = %3.2f Peak = %3.2f Split=%3.2f Base=%3.2f' %
              (coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
    plt.savefig(quicklook_output_dir + '/ch7/' + str(i) + '.png', dpi=300)
    plt.close()

plt.plot(my_binary['amplitude_ch0'].values)
plt.show()