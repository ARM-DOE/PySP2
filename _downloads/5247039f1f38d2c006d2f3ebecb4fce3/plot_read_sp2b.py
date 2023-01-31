"""
Example for plotting a wave in a .sp2b file
-------------------------------------------------------------

"""
import pysp2
import matplotlib.pyplot as plt

my_sp2 = pysp2.io.read_sp2(pysp2.testing.EXAMPLE_SP2B)
my_config = pysp2.io.read_config(pysp2.testing.EXAMPLE_INI)
my_sp2 = pysp2.util.gaussian_fit(my_sp2, my_config)
print(my_sp2)
pysp2.vis.plot_wave(my_sp2, 1100, 0)
plt.show()
