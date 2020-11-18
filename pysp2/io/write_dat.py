""" Function that will write intermediate .dat file for Igor processing"""

import xarray as xr
import pandas as pd
import numpy as np

def write_dat(ds, file_name):
    """
    This writes the .dat files that generate the intermediate parameters used
    by the Igor processing.

    Parameters
    ----------
    ds: xarray Dataset
        The dataset containing the processed signal data.
    file_name: str
        The name of the file to save to.
    """

    ### Convert the dataset to a Pandas DataFrame
    index_label = ["TimeWave", "DateTimeWaveUTC",
                   "DateTimeWave", "EventIndex",
                   "Flag", "Base_ch0", "FtAmp_ch0", "FtPos_ch0",
                   "PkHt_ch0", "PkPos_ch0", "PkFWHM_ch0",
                   "GaussChiSq_ch0", "GaussErrorCode_ch0",
                   "PkStart_ch0", "Base_ch1", "PkHt_ch1",
                   "PkPos_ch1", "PkStart_ch1", "PkEnd_ch1",
                   "PkHalfRise_ch1", "PkHalfDecay_ch1",
                   "Peak2area_ch1", "Base_ch2",
                   "PkHt_ch2", "PkPos_ch2", "PkStart_ch2",
                   "PkEnd_ch2", "PkHalfRise_ch2", "PkHalfDecay_ch2",
                   "Peak2area_ch2", "Base_ch3", "PkHt_ch3",
                   "PkPos_ch3", "PkSplitPos_ch3", "Base_ch4",
                   "FtAmp_ch4", "FtPos_ch4", "PkHt_ch4",
                   "PkPos_ch4", "PkFWHM_ch4", "GaussChiSq_ch4",
                   "GaussErrorCode_ch4", "PkStart_ch4",
                   "Base_ch5", "PkHt_ch5", "PkPos_ch5",
                   "PkStart_ch5", "PkEnd_ch5", "PkHalfRise_ch5",
                   "PkHalfDecay_ch5", "Peak2area_ch5", "Base_ch6",
                   "PkHt_ch6", "PkPos_ch6", "PkStart_ch6", "PkEnd_ch6",
                   "PkHalfRise_ch6", "PkHalfDecay_ch6", "Peak2area_ch6",
                   "Base_ch7", "PkHt_ch7", "PkPos_ch7", "PkSplitPos_ch7",
                   "IncanRatioch5ch6", "IncanPkOffsetch5ch6",
                   "IncanRatioch1ch2", "IncanPkOffsetch1ch2",
                   "ScatRejectKey", "IncanRejectKey"]
    drop_list = []
    for varname in ds.variables.keys():
        if varname not in index_label:
            drop_list.append(varname)

    smaller_ds = ds.drop(drop_list)
    pandas_ds = smaller_ds.to_dataframe()
    sp2_header = ["Instrument Type=SP2\n", "****\n"]

    # Round some entries to fewer figures
    pandas_ds['TimeWave'] = pandas_ds['TimeWave'].map(lambda x: "%.16g" % x)
    pandas_ds['DateTimeWave'] = pandas_ds['DateTimeWave'].map(lambda x: "%.16g" % x)
    pandas_ds['DateTimeWaveUTC'] = pandas_ds['DateTimeWaveUTC'].map(lambda x: "%.16g" % x)
    pandas_ds["GaussChiSq_ch0"] = pandas_ds["GaussChiSq_ch0"].map(lambda x: "%.16g" % x)
    pandas_ds["GaussChiSq_ch4"] = pandas_ds["GaussChiSq_ch4"].map(lambda x: "%.16g" % x)
    pandas_ds["GaussChiSq_ch0"] = pandas_ds["GaussChiSq_ch0"].replace('nan', '')
    pandas_ds["GaussChiSq_ch4"] = pandas_ds["GaussChiSq_ch4"].replace('nan', '')
    pandas_ds["GaussErrorCode_ch0"] = pandas_ds["GaussChiSq_ch0"].map(lambda x: x*0)
    pandas_ds["GaussErrorCode_ch4"] = pandas_ds["GaussChiSq_ch4"].map(lambda x: x*0)
    pandas_ds["IncanRatioch1ch2"] = pandas_ds["IncanRatioch1ch2"].map(lambda x: "%.16g" % x)
    pandas_ds["IncanRatioch5ch6"] = pandas_ds["IncanRatioch5ch6"].map(lambda x: "%.16g" % x)
    pandas_ds["IncanRatioch1ch2"] = pandas_ds["IncanRatioch1ch2"].replace('nan', '')
    pandas_ds["IncanRatioch5ch6"] = pandas_ds["IncanRatioch5ch6"].replace('nan', '')
    with open(file_name, 'w', newline='\n') as f:
        for line_in_header in sp2_header:
            f.write(line_in_header)
        pandas_ds = pandas_ds[index_label]
        #print(pandas_ds)
        pandas_ds.to_csv(f, header=True, index=False, float_format="%.8g", sep='\t', encoding='utf-8')

def write_dat_psds(ds, file_name):
    """
    This writes the .dat files that generate the intermediate parameters used
    by the Igor processing.

    Parameters
    ----------
    ds: xarray Dataset
        The dataset containing the processed signal data.
    file_name: str
        The name of the file to save to.
    """

    pandas_df = pd.DataFrame()
    time = ds.time.values

