"""
Python module for reading SP2-XR Particle-by-Particle (PbP) files.
"""

import os
import re
import warnings

import numpy as np
import pandas as pd

from .read_hk import read_sp2xr_hk_file, _attach_sp2xr_calibration_attrs


# SP2-XR PbP column names mapped to PySP2-style equivalents. Variables on the
# right are intended to look like the output of pysp2.util.gaussian_fit() so
# that downstream calibration utilities can be used with SP2-XR data.
_PBP_RENAME_MAP = {
    'Scatter relPeak': 'PkHt_ch0',
    'Scatter Peak Time': 'PkPos_ch0',
    'Scatter FWHM': 'PkFWHM_ch0',
    'Scatter Transit Time': 'ScatTransitTime',
    'Incand relPeak': 'PkHt_ch1',
    'Incand Peak Time': 'PkPos_ch1',
    'Incand FWHM': 'PkFWHM_ch1',
    'Incand Transit Time': 'IncTransitTime',
    'Incand Delay': 'IncanDelay',
    'Particle Flags': 'ParticleFlags',
    'Particle Time Stamp': 'ParticleTimeStamp',
    'Packet Time Stamp': 'PacketTimeStamp',
    'Time (sec)': 'TimeWave',
    'Dropped Records': 'DroppedRecords',
    'Record Count': 'RecordCount',
    'Record Size': 'RecordSize',
}


def _find_matching_hk_file(pbp_path):
    """
    Locate the SP2-XR HK file that corresponds to a PbP file.

    There is one HK file per acquisition session (always named ``_x0001``),
    while a session may produce many PbP files (``_x0001``, ``_x0002``, ...).
    Both .csv and .zip variants are checked.

    Returns None if no candidate exists on disk.
    """
    base = re.sub(r'PbP', 'hk', pbp_path)
    base = re.sub(r'(_x)\d{4}', r'\g<1>0001', base)
    base = re.sub(r'\.(csv|zip)$', '', base)

    for ext in ('.csv', '.zip'):
        candidate = base + ext
        if os.path.exists(candidate):
            return candidate
    return None


def read_sp2xr_pbp(file_name, keep_firmware_calibration=False):
    """
    Read a SP2-XR Particle-by-Particle (PbP) CSV or ZIP file into an
    xarray Dataset.

    The SP2-XR firmware performs peak fitting on-board, so the PbP file
    already contains per-particle peak heights, positions, FWHMs, transit
    times, and pre-calibrated optical diameter / BC mass. This reader
    maps the peak-fit outputs to PySP2-style names (PkHt_ch0, PkHt_ch1,
    ...) so that the dataset can be fed into PySP2 calibration utilities
    (e.g. pysp2.util.calc_diams_masses) with appropriate SP2-XR
    coefficients.

    By default, the pre-calibrated 'Scatter Size (nm)' and 'Incand Mass
    (fg)' columns are dropped to encourage recalibration with consistent
    curves. Pass keep_firmware_calibration=True to retain them as
    'firmware_ScatterSize_nm' and 'firmware_IncandMass_fg'; in that case
    the scattering / incandescence calibration CSVs ('*_Scatt_*.csv' and
    '*_Incan_*.csv') are auto-located alongside the PbP file and attached
    as dataset attributes for provenance.

    Absolute particle datetimes are derived from the matching SP2-XR HK
    file, which is auto-located in the same directory. Each acquisition
    session writes one HK file (always ``_x0001``) and one or more PbP
    files; the HK provides the wall-clock reference that converts the
    PbP instrument seconds to UTC. If no matching HK file is found, a
    warning is emitted and the dataset is returned without an absolute
    'time' coordinate.

    Parameters
    ----------
    file_name: str
        Path to the PbP CSV or ZIP file.
    keep_firmware_calibration: bool
        If True, retain the SP2-XR firmware's 'Scatter Size (nm)' and
        'Incand Mass (fg)' as 'firmware_ScatterSize_nm' and
        'firmware_IncandMass_fg'. Default False.

    Returns
    -------
    pbp_ds: xarray.Dataset
        Per-particle data indexed by 'event_index'. Variable names follow
        PySP2 conventions for peak fit outputs (PkHt_ch0, PkHt_ch1,
        PkFWHM_ch0, ...), with SP2-XR-specific columns kept under
        descriptive names (ScatTransitTime, IncTransitTime, IncanDelay,
        ParticleFlags).
    """

    # pandas handles both .csv and .zip transparently
    df = pd.read_csv(file_name)

    # Drop the 'Reserved' column (always empty) and, by default, the
    # firmware-calibrated diameter / mass.
    drop_cols = ['Reserved']
    if not keep_firmware_calibration:
        drop_cols += ['Scatter Size (nm)', 'Incand Mass (fg)']
    df = df.drop(columns=[c for c in drop_cols if c in df.columns])

    rename_map = dict(_PBP_RENAME_MAP)
    if keep_firmware_calibration:
        rename_map['Scatter Size (nm)'] = 'firmware_ScatterSize_nm'
        rename_map['Incand Mass (fg)'] = 'firmware_IncandMass_fg'
    # The packet-level 'Flag' shares its name with the per-particle flag
    # column; rename it explicitly to avoid the collision.
    if 'Flag' in df.columns:
        rename_map['Flag'] = 'PacketFlag'
    df = df.rename(columns=rename_map)

    # Use a monotonic event_index as the primary dimension, consistent
    # with read_sp2()/read_sp2xr().
    df.index = pd.Index(np.arange(len(df), dtype='float32'), name='event_index')
    ds = df.to_xarray()
    ds['EventIndex'] = ds['event_index']

    # Attach calibration curve provenance only when the firmware-calibrated
    # size / mass columns are retained.
    if keep_firmware_calibration:
        _attach_sp2xr_calibration_attrs(ds, file_name)

    # Auto-locate the matching HK file to derive absolute datetimes.
    # The filename itself cannot be trusted: it always encodes the session
    # start time, so PbP files past _x0001 would be off by the duration
    # of the preceding files. The HK file provides the only safe anchor.
    hk_file_name = _find_matching_hk_file(file_name)
    if hk_file_name is None:
        warnings.warn(
            "No matching SP2-XR HK file found next to %r; "
            "returning relative instrument seconds only. "
            "Place the HK file alongside the PbP file (with '_x0001') to "
            "get absolute particle datetimes." % file_name,
            UserWarning,
            stacklevel=2,
        )
        return ds

    hk_ds = read_sp2xr_hk_file(hk_file_name)
    t0 = pd.Timestamp(hk_ds['time'].values[0])
    # After read_sp2xr_hk_file strips the unit suffix, the HK 'Time (sec)'
    # column is exposed as 'Time'.
    hk_first_sec = float(hk_ds['Time'].values[0])
    delta_sec = ds['ParticleTimeStamp'].values - hk_first_sec
    time_values = t0 + pd.to_timedelta(delta_sec, unit='s')
    ds = ds.assign_coords(time=('event_index', time_values))

    return ds
