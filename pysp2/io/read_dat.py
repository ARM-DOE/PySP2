""" I/O function for reading .dat files from original IGOR processing."""

import xarray as xr
import pandas as pd
import numpy as np
import act
from glob import glob

def read_dat(file_name, type):
    """
    This reads the .dat files that generate the intermediate parameters used
    by the Igor processing. Wildcards are supported.

    Parameters
    ----------
    file_name: str
        The name of the file to save to. Use a wildcard to open multiple files at once.
    type: str
        This parameter must be one of:
            'particle': Load individual particle timeseries from .dat file
            'conc': Load timeseries of concentrations.
    Returns
    -------
    ds: xarray Dataset
        The xarray dataset to store the parameters in.
    """

    if type.lower() not in ['particle', 'conc']:
        raise ValueError("Invalid input for type, must be either 'particle' or 'conc'!")

    fname = glob(file_name)
    ds_list = []
    for f in fname:
        ds_list.append(act.io.csvfiles.read_csv(f, sep="\t", skiprows=2))

    if type.lower() == 'particle':
        return xr.concat(ds_list, dim='index').sortby('DateTimeWave')
    elif type.lower() == 'conc':
        return xr.concat(ds_list, dim='index').sortby('Start DateTime')




