""" Function that will write intermediate .dat file for Igor processing"""

import pandas as pd
import datetime
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
        pandas_ds.to_csv(f, header=True, index=False, float_format="%.8g", sep='\t', encoding='utf-8')


def write_dat_concs(ds, file_name):
    """
    This writes the .dat files for the mass and number concentrations
    by the Igor processing.

    Parameters
    ----------
    ds: xarray Dataset
        The dataset containing the processed signal data.
    file_name: str
        The name of the file to save to.
    """

    pandas_df = pd.DataFrame()
    cur_date = datetime.datetime(ds.time.dt.year[0], ds.time.dt.month[0], ds.time.dt.day[0])
    try:
        interval = (ds.time.values[1] - ds.time.values[0]) / np.timedelta64(1, 's')
    except IndexError:
        interval = 10
    start_time = ds.time.dt.hour.values * 3600 + ds.time.dt.minute.values * 60 + ds.time.dt.second.values
    end_time = ds.time.dt.hour.values * 3600 + ds.time.dt.minute.values * 60 + \
               ds.time.dt.second.values + interval
    time_wave = (ds.time.values - np.datetime64('1904-01-01T00:00:00')) / np.timedelta64(1, 's')
    pandas_df['Start time'] = start_time
    pandas_df['End time'] = end_time
    pandas_df['Start DateTime'] = time_wave
    pandas_df['End DateTime'] = time_wave + interval
    pandas_df['Scattering conc (#/cm3-STP)'] = ds.NumConcScat.values
    pandas_df['Incandescent conc (#/cm3-STP)'] = ds.NumConcIncan.values
    pandas_df['Total scattering mass conc (ng/m3-STP)'] = ds.MassScat2.values
    pandas_df['Scattering mass conc (ng/m3-STP)'] = ds.MassScat2.values
    pandas_df['Incandescent mass conc (ng/m3-STP)'] = ds.MassIncand2total.values

    pandas_df['Start time'] = pandas_df['Start time'].map(lambda x: '%5d' % x)
    pandas_df['End time'] = pandas_df['End time'].map(lambda x: '%5d' % x)
    pandas_df['Start DateTime'] = pandas_df['Start DateTime'].map(lambda x: '%11d' % x)
    pandas_df['End DateTime'] = pandas_df['End DateTime'].map(lambda x: '%11d' % x)
    pandas_df['Scattering conc (#/cm3-STP)'] = pandas_df['Scattering conc (#/cm3-STP)'].map(lambda x: "%.12g" % x)
    pandas_df['Scattering mass conc (ng/m3-STP)'] = pandas_df['Scattering mass conc (ng/m3-STP)'].map(lambda x: "%.12g" % x)
    pandas_df['Incandescent mass conc (ng/m3-STP)'] = pandas_df['Incandescent mass conc (ng/m3-STP)'].map(
        lambda x: "%.12g" % x)
    pandas_df['Incandescent conc (#/cm3-STP)'] = pandas_df['Incandescent conc (#/cm3-STP)'].map(lambda x: "%.12g" % x)
    pandas_df['Total scattering mass conc (ng/m3-STP)'] = pandas_df['Total scattering mass conc (ng/m3-STP)'].map(lambda x: "%.12g" % x)
    with open(file_name, 'w', newline='\n') as f:
        pandas_df.to_csv(f, header=True, index=False, float_format="%.8g", sep='\t', encoding='utf-8')


def write_dat_concs_arm(ds, file_name, location, lat_lon_string, deltaSize=0.005):
    """
    Write the SP2 size distribution and mass concentration in the
    standard ARM convention.

    Parameters
    ----------
    ds: xarray Dataset
        The dataset with the process timeseries data.
    file_name: str
        The output filename.
    location: str
        The output location.
    deltaSize: float
        The bin width in micrometers.
    """
    my_file = open(file_name, 'w')
    my_file.write("rBC (refractory Black Carbon) concentration at STP (ng/m3)" + 
                  " and number size distribution time series\n")
    my_file.write("Time resolution: 60 sec\n")
    my_file.write("Uncertainty: ~30%\n")
    my_file.write("SP2 Unit: 25")
    my_file.write("rBC mass concentration is reported at STP conditions")
    my_file.write("Location: " + location + '\n')
    my_file.write("Location: " + lat_lon_string + '\n')
    my_file.write("Elevation ~ 300 m asl.\n")
    my_file.write("Calibration: AquaDag with 1.3 factor\n")
    my_file.write("\n\n\n\n\nSP2_dateTime: UTC\n")
    my_file.write("Mass Equivalent Diameters [MED] used for size " + 
                  "distribution (SP2_min; SP2_geo; and SP2_max) are " + 
                  "in units of micrometers\n\n")
    my_file.write("Column naming convention: 'SP2_cnts_X' are the number of" + 
                  " particles in bin number _X. , where _X is the row number\n")
    my_file.write("within the 'SP2_geo' size bin column that contains the mass " + 
                  "equivalent diameter (e.g., SP2_cnts_0 = 0.01 micrometers; " +
                  "SP2_cnts_10 = 0.060 micrometers, etc…).\n")
    my_file.write("The dN/dlogDp data is time-resolved where a given row is " +
                  " associated with the timestamp for that row.\n")
    my_file.write("\nrBC concentration is in units of ng/m^3. " +
                  " Note that the rBC column length is one field.\n")
    my_file.write("shorter than the SP2_datetime column.\n " +
                  " Last time field is not relevant to the rBC time series.\n")
    my_file.write("(see comment below on length of SP2_datetime column)\n")
    my_file.write("\nLengths for SP2_max; SP2_min; SP2_geo are one field longer" + 
                  " then the number of\n")
    my_file.write("SP2_cnts_XX columns .  This is to provide bounds " + 
                  " for image plots (if desired).\n")
    my_file.write("\nLength for SP2_datetime is  one field longer than " + 
                  " that length of the SP2_cnts_XX columns\n")
    my_file.write("This is to provide bounds for image plots (if desired)\n")
    my_file.write("\n\nComments and/or requests are to be directed to Art " +
                  " Sedlacek (sedlacek@bnl.gov)\n")
    out_df = {}

    def time_round(x):
        return "%12.2f" % x
    out_df["SP2_datetime_in_sec"] = ds.TimeWave.values
    out_df["SP2_date"] = ds.time.dt.strftime("%Y/%m/%d")
    out_df["SP2_time"] = ds.time.dt.strftime("%H:%M:%S")
    out_df["SP2_rBC_conc"] = ds.MassIncand2total.values
    num_times = ds.ScatNumEnsembleBC.values.shape[0]
    num_bins = ds.ScatNumEnsembleBC.values.shape[1]
    SpecSizeBins = 0.01 + np.arange(0, num_bins + 1, 1) * deltaSize
    dlogDp = np.diff(np.log10(SpecSizeBins))

    def string_round(x):
        return "%g" % x
    out_df["SP2_Dmin"] = SpecSizeBins[:-1]
    out_df["SP2_Dgeo"] = np.sqrt(SpecSizeBins[:-1] * SpecSizeBins[1:])
    out_df["SP2_Dmax"] = SpecSizeBins[1:]
    for i in range(num_bins):
        out_df["SP2_cnts_%d" % i] = ds.IncanNumEnsemble.values[:, i] / dlogDp[i]
    
    if num_times > num_bins:
        out_df["SP2_Dmin"] = np.pad(
            out_df["SP2_Dmin"], (0, num_times - num_bins), mode='constant', constant_values=np.nan)
        out_df["SP2_Dgeo"] = np.pad(
            out_df["SP2_Dgeo"], (0, num_times - num_bins), mode='constant', constant_values=np.nan)
        out_df["SP2_Dmax"] = np.pad(
            out_df["SP2_Dmax"], (0, num_times - num_bins), mode='constant', constant_values=np.nan)
    elif num_times < num_bins:
        out_df["SP2_date"] = np.pad(
            out_df["SP2_date"], (0, num_bins - num_times), mode='constant', constant_values="")
        out_df["SP2_time"] = np.pad(
            out_df["SP2_time"], (0, num_bins - num_times), mode='constant', constant_values="")
        out_df["SP2_rBC_conc"] = np.pad(
            out_df["SP2_rBC_conc"], (0, num_bins - num_times), mode='constant', constant_values=np.nan)
        for i in range(num_bins):
            out_df["SP2_cnts_%d" % i] = np.pad(
                out_df["SP2_cnts_%d" % i], (0, num_bins - num_times), mode='constant', constant_values=np.nan)
    out_df = pd.DataFrame(out_df)
    out_df["SP2_Dmin"] = out_df["SP2_Dmin"].apply(string_round)
    out_df["SP2_Dgeo"] = out_df["SP2_Dgeo"].apply(string_round)
    out_df["SP2_Dmax"] = out_df["SP2_Dmax"].apply(string_round)
    out_df["SP2_datetime_in_sec"] = out_df["SP2_datetime_in_sec"].apply(time_round)
    out_df = out_df.fillna('')
    out_df.to_csv(my_file, header=True, index=False, float_format="%g", sep='\t', encoding='utf-8')
    my_file.close()
