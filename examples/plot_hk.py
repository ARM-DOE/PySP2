"""
Example for plotting number of particles in housekeeping data
-------------------------------------------------------------

"""
import pysp2
import matplotlib.pyplot as plt

my_hk = pysp2.io.read_hk_file(pysp2.testing.EXAMPLE_HK)
print(my_hk)
# my_hk is a standard xarray dataset, so one can use xarray's plotting
# capabilities in order to make the plot
my_hk['Num of Particles'].plot()
plt.show()
