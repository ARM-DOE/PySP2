"""
pysp2.io
--------

This module contains all of the procedures for reading and writing SP2 data.

.. autosummary::
    :toctree: generated/

    read_hk_file
    get_hk_variable_names
    read_sp2
    read_config
    write_dat
    write_dat_concs
    write_dat_concs_arm
    read_dat
    read_arm_dat
    read_calibration

"""

from .read_hk import read_hk_file, get_hk_variable_names
from .read_hk import read_multi_hk_file
from .read_sp2b import read_sp2
from .read_ini import read_config
from .write_dat import write_dat, write_dat_concs, write_dat_concs_arm
from .read_dat import read_dat, read_arm_dat, read_calibration 
