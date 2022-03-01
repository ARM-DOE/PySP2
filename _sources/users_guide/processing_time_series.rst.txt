Processing timeseries data
==========================

In order to process the particle data into mass and size distributions, the
data must be both filtered of artifacts that can be caused by particle
coincidence in the beam as well as noise. In order to do this, the particles
are filtered out using criteria related to the particle waveform properties.
These settings are a part of the :func:`pysp2.util.DMTGlobals` class.

Using the DMTGlobals class, one can customize the criteria used to
filter out artifacts in the data. The :func:`pysp2.util.DMTGlobals` structure
is initialized with a default set of values with criteria that are commonly
used for processing. The filtering criteria that you can specify are listed
in the table below.

.. list-table:: List of filtering criteria in DMTGlobals
    :header-rows: 1

    * - Variable name
      - Variable description
      - Default value
    * - ScatMaxPeakHt1
      - Maximum peak height for scattering channel 0
      - 60000
    * - ScatMaxPeakHt2
      - Maximum peak height for scattering channel 4
      - 60000
    * - ScatMinPeakHt1
      - Minimum peak height for scattering channel 0
      - 250
    * - ScatMinPeakHt2
      - Minimum peak height for scattering channel 4
      - 250
    * - ScatMinWidth
      - Minimum scattering width
      - 10
    * - ScatMaxWidth
      - Maximum scattering width
      - 90
    * - ScatMinPeakPos
      - Scattering minimum peak position
      - 20
    * - ScatMaxPeakPos
      - Scattering maximum peak position
      - 90
    * - IncanMinPeakHt1
      - Incandescence minimum peak position ch 1
      - 200
    * - IncanMinPeakHt1
      - Incandescence minimum peak position ch 5
      - 200
    * - IncanMaxPeakHt1
      - Incandescence minimum peak position ch 1
      - 60000
    * - IncanMaxPeakHt1
      - Incandescence minimum peak position ch 5
      - 60000
    * - IncanMinWidth
      - Incandescence minimum peak width
      - 5
    * - IncanMaxWidth
      - Incandescence maximum peak with
      - np.inf
    * - IncanMinPeakPos
      - Incandescence minimum peak position
      - 20
    * - IncanMaxPeakPos
      - Incandescence maximum peak position
      - 90
    * - IncanMinPeakRatio
      - The minimum peak ratio between ch5 and 6 and ch1 and 2
      - 0.1
    * - IncanMaxPeakRatio
      - The minimum peak ratio between ch5 and 6 and ch1 and 2
      - 25
    * - IncanMaxPeakOffset
      - The maximum offset between the peaks of ch5/6 or ch1/2
      - 11

In addition, in order to provide the most accurate estimate of particle
masses and sizes, a calibration is typically performed on the SP2 prior to
a field experiment. These calibrations will output coefficients used to
calculate particle sizes and masses in a .cal file. In order to use the
specified calibration in the mass calculations, you simply enter the
.cal file as the parameter to DMTGlobals when initializing your
structure:

.. code-block:: python

    global_settings = pysp2.util.DMTGlobals('my_calibration.cal')

The :code:`global_settings` are then used as an input to the particle size
and mass calculations:

.. code-block:: python

    particles = pysp2.util.calc_diams_masses(my_ds, globals=globals)

After the particle sizes are calculated, you need to combine both the
housekeeping and configuration data with these particle masses in order to
calculate the particle size distributions. An example snippet that
demonstrates the entire process is below:

.. code-block:: python

    def test_psds():
        my_sp2b = pysp2.io.read_sp2(pysp2.testing.EXAMPLE_SP2B)
        my_ini = pysp2.io.read_config(pysp2.testing.EXAMPLE_INI)
        my_binary = pysp2.util.gaussian_fit(my_sp2b, my_ini)
        my_hk = pysp2.io.read_hk_file(pysp2.testing.EXAMPLE_HK)
        my_binary = pysp2.util.calc_diams_masses(my_binary)
        my_psds = pysp2.util.process_psds(my_binary, my_hk, my_ini)

The my_psds structure contains an xarray dataset with the following
timeseries values:

.. list-table:: Particle number and mass size distributions
   :header-rows: 1

   * - time
     - Time in seconds since the epoch
   * - TimeWave
     - Time in seconds since Jan 1, 1904 at midnight
   * - NumConcIncan
     - Number concentration of incandescence particles [cm-3]
   * - NumConcIncanScat
     - Number concentration of incandescence particles also detected by scattering channel [cm-3]
   * - NumConcTotal
     - Total number concentration of all particles [cm-3]
   * - NumConcScatSat
     - Number concentration of scattering particles that saturate the channel [cm-3]
   * - NumConcScat
     - Number concentration of scattering particles excluding saturated signals [cm-3]
   * - MassIncand2
     - Mass concentration of incandescence particles [ng cm-3]
   * - MassIncand2Sat
     - Mass concentration of particles detected by both incandescence and scattering channels [nm cm-3]
   * - MassIncand2total
     - Total mass concentration of particles [ng cm-3]
   * - MassScat2
     - Mass concentration of particles detected by scattering channel
   * - ScatNumEnsemble
     - Particle size distribution of scattering particles [cm-3 per bin]
   * - ScatMassEnsemble
     - Particle mass size distribution of scattering particles [ng cm-3 per bin]
   * - IncanNumEnsemble
     - Incandescence particle size distribution [cm-3 per bin]
   * - IncanMassEnsemble
     - Incandescence particle mass distribution [ng cm-3 per bin]
   * - ScatNumEnsembleBC
     - Scattering particle size distribution assuming black carbon [ng cm-3 per bin]
   * - NumFracBC
     - Number fraction of black carbon particles

These data can then be either save to netCDF by using the .to_netcdf()
directive from xarray, or to a .dat file using :func:`pysp2.io.write_dat_concs`.

PySP2 also supports the loading of .dat files generated by the IGOR code
using the :func:`pysp2.io.read_dat` function.



