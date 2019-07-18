"""
Python module for reading housekeeping files
"""

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
    hk_df: pandas DataFrame
        The housekeeping information in a pandas DataFrame
    """

    my_df = pd.read_csv(file_name, delimiter="\t")
    my_df = my_df.set_index('Time (sec)')
    return my_df

def get_hk_variable_names(my_df):
    """
    This procedure will return al ist of variables in the
    housekeeping file.

    Parameters
    ----------
    my_df: pandas DataFrame
        The dataframe to get the variable names from

    Returns
    -------
    var_names: list
        The names of each variable in the file.
    """
    return [my_str for my_str in my_df.keys()]

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
    my_df: pandas dataframe
        The pandas DataFrame containing the data loaded.
    """

    the_list = []
    file_list = glob.glob(file_path)

    for f in file_list:
        df = read_hk_file(f)
        the_list.append(df)

    return pd.concat(df, axis=0)