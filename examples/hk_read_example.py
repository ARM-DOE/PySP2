import pysp2
import matplotlib.pyplot as plt

my_hk_file = pysp2.io.read_hk_file('/home/rjackson/Downloads/SP2Data4ANL/20190518000000.hk')
print(pysp2.io.get_hk_variable_names(my_hk_file))
pysp2.vis.plot_hk_timeseries(my_hk_file, ['YAG Power (V)',
                                          'Sample Flow LFE (vccm)'],
                             linewidth=2)
plt.show()