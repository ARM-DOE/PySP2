# PySP2

This is a python package for processing and visualizing SP2 data. It is based off of the IGOR code developed by Droplet Measurement Technologies. 
PySP2 currently supports processing all raw SP2 data (.sp2b, .hk, .ini) files into mass and number size distributions. 
It can plot individual waveforms as well as particle size distributions. Its file format is based off of the standard
provided by the [Atmospheric Community Toolkit](https://arm-doe.github.io/ACT) which is based around [xarray](https://xarray.pydata.org). 
PySP2 is currently used by the Department of Energy Atmospheric Radiation Measurment (ARM) Facility in order to process SP2 observations during field deployments such
as TRACER and SAIL.

[Plot of a waveform](https://arm-doe.github.io/PySP2/_images/sphx_glr_plot_read_sp2b_001.png "Plot of a waveform")
# Important links

Documentation: https://arm-doe.github.io/PySP2/

Examples: https://arm-doe.github.io/PySP2/source/auto_examples/plot_read_sp2b.html

# Installation

PySP2 is pip

# References

For more information about how particle sizes are derived from SP2 signals using calibration, please consult the following references:

M. Irwin, Y. Kondo, N. Moteki & T. Miyakawa: Evaluation of a Heated-Inlet for Calibration of the SP2, Aerosol Science and Technology, 47:8, 895-905, DOI: 10.1080/02786826.2013.800187, 2013

Gysel, M., Laborde, M., Olfert, J. S., Subramanian, R., and Gröhn, A. J.: Effective density of Aquadag and fullerene soot black carbon reference materials used for SP2 calibration, Atmos. Meas. Tech., 4, 2851–2858, https://doi.org/10.5194/amt-4-2851-2011, 2011. 
