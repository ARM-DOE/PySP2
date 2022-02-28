============
User's Guide
============

**Installing**

In order to install PySP2, one can use either pip or anaconda. In order
to pip install pysp2, simply execute:

::

  pip install pysp2


For anaconda users, the command is:

::

   conda install -c conda-forge pysp2


**Usage**

PySP2 depends on `xarray <https://xarray.pydata.org>`_ and the
`Atmospheric Community Toolkit <https://ARM-DOE.github.io/ACT>`_ in order
to handle its data structures. All data structures used by PySP2 are in
standard `xarray <https://xarray.pydata.org>`_ format. Therefore, for
examples on how to analyze `xarray <https://xarray.pydata.org>`_ datasets,
we recommend consulting the `xarray User Guide <https://docs.xarray.dev/en/stable/user-guide/index.html>`_.

The PySP2 user guide shows basic examples of how to use PySP2 to process raw .sp2b files,
view particle waveforms, housekeeping data, as well as for processing particle
size distribution timeseries data.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   processing_raw_file
   particle_information
   view_hk_data
   processing_time_series


