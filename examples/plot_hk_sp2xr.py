"""
Example for plotting selected variables from a SP2-XR housekeeping file
-----------------------------------------------------------------------

"""
import pysp2
import matplotlib.pyplot as plt

my_hk = pysp2.io.read_sp2xr_hk_file(pysp2.testing.EXAMPLE_SP2XR_HK)
print(my_hk)

# A few representative housekeeping channels: temperature, flow, pressure,
# and particle event count
variables = [
    'Laser TEC Temp',
    'Sample Flow Controller Read',
    'Cavity Pressure',
    'Threshold Crossing Events',
]

fig, axes = plt.subplots(len(variables), 1, sharex=True, figsize=(8, 8))
for ax, var in zip(axes, variables):
    my_hk[var].plot(ax=ax)
    ax.set_title(var)
plt.tight_layout()
plt.show()
