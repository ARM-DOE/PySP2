"""
Python module for reading housekeeping files
"""

import xarray as xr
import act
import datetime
import os
import re
import warnings
import numpy as np
import pandas as pd

from glob import glob

# Regex for SP2-XR firmware-computed PSD bin columns. The number of bins per
# channel is not assumed; it is determined from however many columns match
# the pattern in the file at read time.
_SP2XR_BIN_COL_REGEX = re.compile(r'^(?P<channel>Scatter|Incand) Bin (?P<idx>\d+)$')

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

    my_df = act.io.read_csv(file_name, sep="\t")

    # Parse time from filename
    start_time = pd.Timestamp('1904-01-01')
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


def read_sp2xr_hk_file(file_name):
    """
    This procedure will read in an SP2-XR housekeeping file (.csv or .zip)
    and store the timeseries data into an xarray Dataset.

    The SP2-XR housekeeping file uses comma-separated values (rather than
    tab-separated like the original SP2). The 'Time Stamp (UTC sec)' column
    is used to build the time coordinate, sharing the same 1904-01-01 epoch
    as the original SP2 'Timestamp' column.

    The SP2-XR firmware-computed PSD bin columns (named like 'Scatter Bin N'
    and 'Incand Bin N') are excluded from the output to keep the variable
    set similar to read_hk_file(). Use read_sp2xr_hk_psd() to access them
    as a 2D dataset.

    The retained 'Scattering Mass Conc' and 'Incand Mass Conc' columns are
    firmware-computed using the on-board calibration curves. To preserve
    provenance, the calibration CSV files alongside the HK file (matching
    '*_Scatt_*.csv' and '*_Incan_*.csv') are auto-located and attached as
    dataset attributes; a warning is emitted for any missing file.

    Parameters
    ----------
    file_name: str
        The file name to read in (.csv or .zip).

    Returns
    -------
    hk_ds: xarray.Dataset
        The housekeeping information as an xarray Dataset indexed by time,
        with the firmware-computed PSD bin columns removed.
    """

    my_df = act.io.read_csv(file_name, sep=",")

    # Parse time from the UTC seconds column (same 1904 epoch as SP2)
    start_time = pd.Timestamp('1904-01-01')
    my_df = my_df.rename({'index': 'time'})
    my_df['time'] = np.array([
        start_time + datetime.timedelta(seconds=x)
        for x in my_df['Time Stamp (UTC sec)'].values
    ])
    my_df['time'].attrs['units'] = "datetime"
    my_df['time'].attrs['long_name'] = "Time [SP2-XR time]"

    # Drop firmware-computed PSD bin columns
    bin_vars = [v for v in my_df.data_vars
                if _SP2XR_BIN_COL_REGEX.match(str(v))]
    if bin_vars:
        my_df = my_df.drop_vars(bin_vars)

    # Extract units from column names following the "Name (units)" convention
    for vars in my_df.variables.keys():
        splits = vars.split("(")
        try:
            units = splits[1][:-1]
            my_df[vars].attrs['units'] = units
            my_df[vars].attrs['long_name'] = vars
            my_df = my_df.rename({vars: splits[0][:-1]})
        except (IndexError, ValueError):
            continue

    # Attach calibration curve provenance (used by the firmware to derive
    # the retained Mass Conc columns).
    _attach_sp2xr_calibration_attrs(my_df, file_name)

    return my_df


def _find_sp2xr_calibration_files(hk_path):
    """
    Locate the SP2-XR scattering and incandescence calibration CSVs that
    sit alongside an HK file. The convention is filenames containing
    '_Scatt_' (scattering) and '_Incan_' (incandescence) in the same
    directory as the HK file.

    Returns
    -------
    (scatt_path, incan_path): tuple
        Each element is the absolute path to the calibration CSV, or None
        if no matching file is present.
    """
    dir_path = os.path.dirname(os.path.abspath(hk_path))
    scatt = sorted(glob(os.path.join(dir_path, '*_Scatt_*.csv')))
    incan = sorted(glob(os.path.join(dir_path, '*_Incan_*.csv')))
    return (scatt[0] if scatt else None,
            incan[0] if incan else None)


def _load_sp2xr_calibration_curve(csv_path):
    """
    Load a SP2-XR calibration CSV (no header, two columns: physical value
    and instrument signal). Empty trailing rows are dropped.
    """
    df = pd.read_csv(csv_path, header=None).dropna()
    return df.iloc[:, 0].to_numpy(), df.iloc[:, 1].to_numpy()


def _attach_sp2xr_calibration_attrs(ds, file_path):
    """
    Auto-locate scattering and incandescence calibration CSVs alongside
    file_path and attach them as dataset attributes. Emits a UserWarning
    for each calibration file that cannot be found.

    Attached attributes (when available):
        ScatCalibrationFile, ScatCalibration_Diameter_nm,
        ScatCalibration_Signal,
        IncanCalibrationFile, IncanCalibration_Mass_fg,
        IncanCalibration_Signal.
    """
    scatt_cal, incan_cal = _find_sp2xr_calibration_files(file_path)

    if scatt_cal is not None:
        diam_nm, signal = _load_sp2xr_calibration_curve(scatt_cal)
        ds.attrs['ScatCalibrationFile'] = scatt_cal
        ds.attrs['ScatCalibration_Diameter_nm'] = diam_nm
        ds.attrs['ScatCalibration_Signal'] = signal
    else:
        warnings.warn(
            "No SP2-XR scattering calibration file ('*_Scatt_*.csv') found "
            "next to %r; ScatCalibration_* attributes will not be attached."
            % file_path,
            UserWarning,
            stacklevel=3,
        )

    if incan_cal is not None:
        mass_fg, signal = _load_sp2xr_calibration_curve(incan_cal)
        ds.attrs['IncanCalibrationFile'] = incan_cal
        ds.attrs['IncanCalibration_Mass_fg'] = mass_fg
        ds.attrs['IncanCalibration_Signal'] = signal
    else:
        warnings.warn(
            "No SP2-XR incandescence calibration file ('*_Incan_*.csv') "
            "found next to %r; IncanCalibration_* attributes will not be "
            "attached." % file_path,
            UserWarning,
            stacklevel=3,
        )


def read_sp2xr_hk_psd(file_name):
    """
    Read the SP2-XR firmware-computed PSD histograms from a HK file.

    The on-board firmware accumulates per-second particle counts into
    size/mass bins ('Scatter Bin N' for scattering size, 'Incand Bin N'
    for incandescence mass) using its own calibration. This procedure
    extracts those columns into a 2D dataset shaped (time, num_bins),
    matching the convention used by pysp2.util.process_psds() for
    PySP2-computed distributions (ScatNumEnsemble, IncanNumEnsemble).

    Calibration curves for scattering diameter and incandescence mass are
    auto-located alongside the HK file (files matching '*_Scatt_*.csv' and
    '*_Incan_*.csv' respectively) and attached as dataset attributes when
    present. A warning is issued for any missing calibration file.

    Parameters
    ----------
    file_name: str
        The HK file name to read in (.csv or .zip).

    Returns
    -------
    psd_ds: xarray.Dataset
        Dataset with dimensions (time, num_bins) and data variables
        ScatNumEnsemble and IncanNumEnsemble (raw counts per bin per
        time). Sample volume must be applied to convert to concentration.
    """

    my_df = act.io.read_csv(file_name, sep=",")

    start_time = pd.Timestamp('1904-01-01')
    my_df = my_df.rename({'index': 'time'})
    my_df['time'] = np.array([
        start_time + datetime.timedelta(seconds=x)
        for x in my_df['Time Stamp (UTC sec)'].values
    ])
    my_df['time'].attrs['units'] = "datetime"
    my_df['time'].attrs['long_name'] = "Time [SP2-XR time]"

    # Discover bin columns dynamically (do not assume a fixed bin count)
    scat_cols, incan_cols = {}, {}
    for var in list(my_df.data_vars):
        m = _SP2XR_BIN_COL_REGEX.match(str(var))
        if m is None:
            continue
        idx = int(m.group('idx'))
        if m.group('channel') == 'Scatter':
            scat_cols[idx] = str(var)
        else:
            incan_cols[idx] = str(var)

    if not scat_cols and not incan_cols:
        raise ValueError(
            "No SP2-XR firmware PSD bin columns ('Scatter Bin N', "
            "'Incand Bin N') found in %r." % file_name)

    def _stack(cols_by_idx):
        ordered = [cols_by_idx[i] for i in sorted(cols_by_idx)]
        return np.stack([my_df[c].values for c in ordered], axis=-1)

    note = ('Raw particle counts per bin per time. To get concentration in '
            'cts/cm^3 divide by the sample volume per record '
            '(Sample Flow Controller Read (vccm) / 60).')

    psd_ds = xr.Dataset(coords={'time': my_df['time']})

    if scat_cols:
        psd_ds['ScatNumEnsemble'] = xr.DataArray(
            _stack(scat_cols), dims=('time', 'num_bins'),
            attrs={
                'long_name':
                    'Scattering number distribution (SP2-XR firmware, raw counts)',
                'standard_name': 'scattering_number_distribution',
                'units': 'count',
                'source': 'SP2-XR firmware on-board calculation',
                'note': note,
            })

    if incan_cols:
        psd_ds['IncanNumEnsemble'] = xr.DataArray(
            _stack(incan_cols), dims=('time', 'num_bins'),
            attrs={
                'long_name':
                    'Incandescence number distribution (SP2-XR firmware, raw counts)',
                'standard_name': 'incandescence_number_distribution',
                'units': 'count',
                'source': 'SP2-XR firmware on-board calculation',
                'note': note,
            })

    psd_ds.attrs['source'] = 'SP2-XR firmware'

    _attach_sp2xr_calibration_attrs(psd_ds, file_name)

    return psd_ds


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
        Examples: 
            Read all .hk files in one directoy:
                my_hk = pysp2.io.read_multi_hk_file('/path/to/directory/*.hk')
            Read all .hk files and check in the subdirectories as well.
                my_hk = pysp2.io.read_multi_hk_file('/path/to/directory/**/*.hk')
                
    Returns
    -------
    my_df: xarray.Dataset
        The xarray Dataset containing the data loaded.
    """

    the_list = []
    file_list = glob(file_path)

    for f in file_list:
        df = read_hk_file(f)
        the_list.append(df)

    return xr.concat(the_list, dim='time').sortby('time')
