"""
Python module for reading housekeeping files
"""

import xarray as xr
import act
import datetime
import os
import numpy as np
import pandas as pd

from glob import glob


def read_hk_file(file_name):
    """
    This procedure will read in an SP2 housekeeping file and then
    store the timeseries data into a pandas DataFrame.

    Parameters
    ----------
    file_name: str
        The file name to read in

    Returns
    -------
    hk_df: pandas.DataFrame
        The housekeeping information in a pandas DataFrame
    """

    my_df = act.io.csvfiles.read_csv(file_name, sep="\t")
    # Parse time from filename
    start_time = pd.Timestamp('1904-01-01')
    my_df = my_df.set_index({'index': 'Time (sec)'})
    my_df = my_df.rename({'index': 'time'})
    my_df['time'] = np.array([start_time + datetime.timedelta(seconds=x) for x in my_df['Timestamp'].values])
    my_df['time'].attrs['units'] = "datetime"
    my_df['time'].attrs['long_name'] = "Time [SP2 time]"
    for vars in my_df.variables.keys():
        splits = vars.split("(")
        try:
            units = splits[1][:-1]
            my_df[vars].attrs['units'] = units
            my_df[vars].attrs['long_name'] = vars
            my_df = my_df.rename({vars: splits[0][:-1]})
        except (IndexError, ValueError):
            continue

    return my_df


def get_hk_variable_names(my_df):
    """
    This procedure will return al ist of variables in the
    housekeeping file.

    Parameters
    ----------
    my_df: xarray.Dataset
        The dataframe to get the variable names from

    Returns
    -------
    var_names: list
        The names of each variable in the file.
    """
    return [my_str for my_str in my_df.variables.keys()]


def read_multi_hk_file(file_path):
    """
    This procedure will read multiple housekeeping files
    and then concatenate them into a single pandas
    DataFrame

    Parameters
    ----------
    file_path: str
        The path (with wildcards) to the housekeeping files.

    Returns
    -------
    my_df: xarray.Dataset
        The xarray Dataset containing the data loaded.
    """

    the_list = []
    file_list = glob.glob(file_path)

    for f in file_list:
        df = read_hk_file(f)
        the_list.append(df)

    return xr.concat(df, dim='index')
