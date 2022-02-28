"""
Example on plotting particle concentration data timeseries
----------------------------------------------------------

"""
import matplotlib.pyplot as plt
import pysp2


my_sp2b = pysp2.io.read_sp2(pysp2.testing.EXAMPLE_SP2B)
my_ini = pysp2.io.read_config(pysp2.testing.EXAMPLE_INI)
my_binary = pysp2.util.gaussian_fit(my_sp2b, my_ini)
my_hk = pysp2.io.read_hk_file(pysp2.testing.EXAMPLE_HK)
my_binary = pysp2.util.calc_diams_masses(my_binary)
my_psds = pysp2.util.process_psds(my_binary, my_hk, my_ini)

my_psds['NumConcIncan'].plot()
plt.show()
