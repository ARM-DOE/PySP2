"""
Example for plotting a waveform in a SP2-XR .sp2b file
------------------------------------------------------

"""
import pysp2
import matplotlib.pyplot as plt

my_sp2xr = pysp2.io.read_sp2xr(pysp2.testing.EXAMPLE_SP2XR_SP2B)
print(my_sp2xr)

# Plot ch0 and ch1 waveforms for the first particle
fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
my_sp2xr['Data_ch0'].isel(event_index=0).plot(ax=ax0)
ax0.set_title('Particle 0 - Data_ch0')
my_sp2xr['Data_ch1'].isel(event_index=0).plot(ax=ax1)
ax1.set_title('Particle 0 - Data_ch1')
plt.tight_layout()
plt.show()
