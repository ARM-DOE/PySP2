import pysp2
import matplotlib.pyplot as plt
import act

my_hk_file = pysp2.io.read_hk_file('C:\\Users\\rjackson\\Downloads\\SP2Data4ANL\\20190518000000.hk')
print(my_hk_file)
Tseries = act.plotting.TimeSeriesDisplay(my_hk_file, subplot_shape=(2,2))
Tseries.plot('YAG Power', subplot_index=(0,0))
Tseries.plot('Incand. Conc.', subplot_index=(1,0))
Tseries.plot('Laser Temp', subplot_index=(0,1))
plt.show()