import numpy as np
import time
import dask.bag as db

from scipy.optimize import curve_fit
from scipy.signal import find_peaks
#from scipy.stats import chisquare
from .DMTGlobals import DMTGlobals


def _do_fit_records(my_ds, i, num_trig_pts, debug=True):
    if debug and i % 100 == 0:
        print("Processing record %d" % i)

    FtAmp = np.zeros(8)
    FtPos = np.zeros(8)
    Base = np.zeros(8)
    PeakHeight = np.zeros(8)
    PeakPos = np.zeros(8)
    GaussChiSq = np.zeros(8)
    PeakStart = np.zeros(8)
    PeakEnd = np.zeros(8)
    Width = np.zeros(8)
    HalfRise = np.zeros(8)
    HalfDecay = np.zeros(8)
    Peak2Area = np.zeros(8)
    Error = np.zeros(8)

    # Do channel 0 first
    coeffs = _fit_record_gaussian(my_ds, i)
    FtAmp[0] = coeffs['amplitude']
    Base[0] = coeffs['base']
    PeakHeight[0] = coeffs['height']
    PeakPos[0] = coeffs['pos']
    FtPos[0] = coeffs['peakpos']
    PeakStart[0] = coeffs['start']
    GaussChiSq[0] = coeffs['chi2']
    Width[0] = coeffs['width']
    Error[0] = coeffs['error']

    for chn in [1, 2, 5, 6]:
        coeffs = _fit_record_incan_ave_base(my_ds, i, chn, num_trig_pts)
        Base[chn] = coeffs['base']
        PeakHeight[chn] = coeffs['height']
        PeakPos[chn] = coeffs['pos']
        PeakStart[chn] = coeffs['start']
        PeakEnd[chn] = coeffs['end']
        HalfRise[chn] = coeffs['half_rise']
        HalfDecay[chn] = coeffs['half_decay']
        Peak2Area[chn] = coeffs['peak2area']

    for chn in [3, 7]:
        coeffs = _split_scatter_fit(my_ds, i, chn)
        Base[chn] = coeffs['base']
        PeakHeight[chn] = coeffs['height']
        PeakPos[chn] = coeffs['pos']
        PeakStart[chn] = coeffs['start']

    # Channel 4 is a special snowflake
    coeffs = _gaussian_sat_fit(my_ds, i)
    FtAmp[4] = coeffs['fitamplitude']
    Base[4] = coeffs['base']
    PeakHeight[4] = coeffs['height']
    PeakPos[4] = coeffs['pos']
    FtPos[4] = coeffs['fitpos']
    PeakStart[4] = coeffs['start']
    GaussChiSq[4] = coeffs['chi2']
    Width[4] = coeffs['width']
    Error[4] = coeffs['error']
    IncanPkOffsetCh1Ch2 = PeakPos[2] - PeakPos[1]
    IncanPkOffsetCh5Ch6 = PeakPos[6] - PeakPos[5]

    # Let's now do the incandescence ratio
    IncanRatioCh1Ch2 = _calc_incan_ratio(my_ds, i, 1, 2, PeakPos, HalfDecay, Base)
    IncanRatioCh5Ch6 = _calc_incan_ratio(my_ds, i, 5, 6, PeakPos, HalfDecay, Base)
    return (FtAmp, FtPos, Base, PeakHeight, PeakPos, GaussChiSq, PeakStart,
            PeakEnd, Width, HalfRise, HalfDecay, Peak2Area, IncanPkOffsetCh1Ch2, IncanPkOffsetCh5Ch6,
            Error, IncanRatioCh1Ch2, IncanRatioCh5Ch6)

def _calc_incan_ratio(my_ds, record_number, ch1, ch2, PeakPos, halfDecay, Base):
    data_ch1 = my_ds['Data_ch' + str(ch1)].isel(event_index=record_number).values
    data_ch2 = my_ds['Data_ch' + str(ch2)].isel(event_index=record_number).values
    ratio = np.nan
    numRatio = 1
    if np.isfinite(PeakPos[ch1]) and np.isfinite(halfDecay[ch1]):
        start = int(PeakPos[ch1])
        stop = int(halfDecay[ch1])+1
        data_ch2 = data_ch2[start:stop]
        data_ch1 = data_ch1[start:stop]
        where_sum = np.where(data_ch2 > Base[ch2])
        ratio = np.sum((data_ch1[where_sum] - Base[ch1])/(data_ch2[where_sum] - Base[ch2]))
        numRatio = len(where_sum[0])

        if numRatio == 0:
            numRatio = 1
            ratio = np.nan
    return ratio/numRatio

def chisquare(obs, f_exp):
    return np.sum((obs - f_exp)**2)

def gaussian_fit(my_ds, config, parallel=False, num_records=None):
    """
    Does Gaussian fitting for each wave in the dataset.
    This will do the fitting for channel 0 only.

    Parameters
    ----------
    my_ds: xarray Dataset
        Raw SP2 binary dataset
    config: ConfigParser object
        The configuration loaded from the INI file.
    parallel: bool
        If true, use dask to enable parallelism
    num_records: int or None
        Only process first num_records datapoints. Set to
        None to process all records.

    Returns
    -------
    wave_ds: xarray Dataset
        Dataset with gaussian fits
    """
    if num_records is None:
        num_records = len(my_ds.Res8.values)

    num_trig_pts = int(config['Acquisition']['Pre-Trig Points'])
    start_time = time.time()
    if not parallel:
        proc_records = []
        for i in range(num_records):
            proc_records.append(_do_fit_records(my_ds, i, num_trig_pts))
    else:
        fit_record = lambda x: _do_fit_records(my_ds, x, num_trig_pts)
        the_bag = db.from_sequence(range(num_records))
        proc_records = the_bag.map(fit_record).compute()

    FtAmp = np.stack([x[0] for x in proc_records])
    FtPos = np.stack([x[1] for x in proc_records])
    Base = np.stack([x[2] for x in proc_records])
    PeakHeight = np.stack([x[3] for x in proc_records])
    PeakPos = np.stack([x[4] for x in proc_records])
    GaussChiSq = np.stack([x[5] for x in proc_records])
    PeakStart = np.stack([x[6] for x in proc_records])
    PeakEnd = np.stack([x[7] for x in proc_records])
    Width = np.stack([x[8] for x in proc_records])
    HalfRise = np.stack([x[9] for x in proc_records])
    HalfDecay = np.stack([x[10] for x in proc_records])
    Peak2Area = np.stack([x[11] for x in proc_records])
    IncanPkOffsetch1ch2 = np.array([x[12] for x in proc_records])
    IncanPkOffsetch5ch6 = np.array([x[13] for x in proc_records])
    Error = np.array([x[14] for x in proc_records])
    IncanRatioch1ch2 = np.array([x[15] for x in proc_records])
    IncanRatioch5ch6 = np.array([x[16] for x in proc_records])

    # Channel 0
    i = 0
    my_ds['FtAmp_ch' + str(i)] = (('event_index'), FtAmp[:,i])
    my_ds['FtAmp_ch' + str(i)].attrs["long_name"] = "Fit Amplitude for channel %d" % i
    my_ds['FtAmp_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['FtPos_ch' + str(i)] = (('event_index'), FtPos[:,i])
    my_ds['FtPos_ch' + str(i)].attrs["long_name"] = "Fit Position for channel %d" % i
    my_ds['FtPos_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['Base_ch' + str(i)] = (('event_index'), Base[:, i])
    my_ds['Base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
    my_ds['Base_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkHt_ch' + str(i)] = (('event_index'), PeakHeight[:, i])
    my_ds['PkHt_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
    my_ds['PkHt_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkFWHM_ch' + str(i)] = (('event_index'), Width[:, i])
    my_ds['PkFWHM_ch' + str(i)].attrs["long_name"] = "Width for channel %d" % i
    my_ds['PkFWHM_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkPos_ch' + str(i)] = (('event_index'), PeakPos[:, i])
    my_ds['PkPos_ch' + str(i)].attrs["long_name"] = "Peak position for channel %d" % i
    my_ds['PkPos_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkStart_ch' + str(i)] = (('event_index'), PeakStart[:, i])
    my_ds['PkStart_ch' + str(i)].attrs["long_name"] = "Peak start for channel %d" % i
    my_ds['PkStart_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['GaussChiSq_ch' + str(i)] = (('event_index'), GaussChiSq[:, i])
    my_ds['GaussChiSq_ch' + str(i)].attrs["long_name"] = "Chisquare value for channel %d" % i
    my_ds['GaussChiSq_ch' + str(i)].attrs["_FillValue"] = np.nan
    
    # Channels 1, 2, 6, 7
    for i in [1, 2, 5, 6]:
        my_ds['FtAmp_ch' + str(i)] = (('event_index'), FtAmp[:, i])
        my_ds['FtAmp_ch' + str(i)].attrs["long_name"] = "Amplitude for channel %d" % i
        my_ds['FtAmp_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['Base_ch' + str(i)] = (('event_index'), Base[:, i])
        my_ds['Base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
        my_ds['Base_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkHt_ch' + str(i)] = (('event_index'), PeakHeight[:, i])
        my_ds['PkHt_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
        my_ds['PkHt_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkHalfRise_ch' + str(i)] = (('event_index'), HalfRise[:, i])
        my_ds['PkHalfRise_ch' + str(i)].attrs["long_name"] = "Point where rise is at 1/2 height for channel %d" % i
        my_ds['PkHalfRise_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['Peak2area_ch' + str(i)] = (('event_index'), Peak2Area[:, i])
        my_ds['Peak2area_ch' + str(i)].attrs["long_name"] = "Peak 2 area for channel %d" % i
        my_ds['Peak2area_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkHalfDecay_ch' + str(i)] = (('event_index'), HalfDecay[:, i])
        my_ds['PkHalfDecay_ch' + str(i)].attrs["long_name"] = "Point where decay is at 1/2 height for channel %d" % i
        my_ds['PkHalfDecay_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkPos_ch' + str(i)] = (('event_index'), PeakPos[:, i])
        my_ds['PkPos_ch' + str(i)].attrs["long_name"] = "Peak position for channel %d" % i
        my_ds['PkPos_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkStart_ch' + str(i)] = (('event_index'), PeakStart[:, i])
        my_ds['PkStart_ch' + str(i)].attrs["long_name"] = "Peak start for channel %d" % i
        my_ds['PkStart_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkEnd_ch' + str(i)] = (('event_index'), PeakEnd[:, i])
        my_ds['PkEnd_ch' + str(i)].attrs["long_name"] = "Peak end for channel %d" % i
        my_ds['PkEnd_ch' + str(i)].attrs["_FillValue"] = np.nan


    for i in [3, 7]:
        my_ds['FtAmp_ch' + str(i)] = (('event_index'), FtAmp[:, i])
        my_ds['FtAmp_ch' + str(i)].attrs["long_name"] = "Amplitude for channel %d" % i
        my_ds['FtAmp_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['Base_ch' + str(i)] = (('event_index'), Base[:, i])
        my_ds['Base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
        my_ds['Base_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkHt_ch' + str(i)] = (('event_index'), PeakHeight[:, i])
        my_ds['PkHt_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
        my_ds['PkHt_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkFWHM_ch' + str(i)] = (('event_index'), Width[:, i])
        my_ds['PkFWHM_ch' + str(i)].attrs["long_name"] = "Width for channel %d" % i
        my_ds['PkFWHM_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkSplitPos_ch' + str(i)] = (('event_index'), PeakStart[:, i])
        my_ds['PkSplitPos_ch' + str(i)].attrs["long_name"] = "Peak start position for channel %d" % i
        my_ds['PkSplitPos_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkPos_ch' + str(i)] = (('event_index'), PeakPos[:, i])
        my_ds['PkPos_ch' + str(i)].attrs["long_name"] = "Peak split position for channel %d" % i
        my_ds['PkPos_ch' + str(i)].attrs["_FillValue"] = np.nan

    # Channel 4
    i = 4
    my_ds['FtAmp_ch' + str(i)] = (('event_index'), FtAmp[:, i])
    my_ds['FtAmp_ch' + str(i)].attrs["long_name"] = "Amplitude for channel %d" % i
    my_ds['FtAmp_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['Base_ch' + str(i)] = (('event_index'), Base[:, i])
    my_ds['Base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
    my_ds['Base_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkHt_ch' + str(i)] = (('event_index'), PeakHeight[:, i])
    my_ds['PkHt_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
    my_ds['PkHt_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkFWHM_ch' + str(i)] = (('event_index'), Width[:, i])
    my_ds['PkFWHM_ch' + str(i)].attrs["long_name"] = "Width for channel %d" % i
    my_ds['PkFWHM_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkPos_ch' + str(i)] = (('event_index'), PeakPos[:, i])
    my_ds['PkPos_ch' + str(i)].attrs["long_name"] = "Peak position for channel %d" % i
    my_ds['PkPos_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['FtPos_ch' + str(i)] = (('event_index'), FtPos[:, i])
    my_ds['FtPos_ch' + str(i)].attrs["long_name"] = "Fit position for channel %d" % i
    my_ds['FtPos_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkStart_ch' + str(i)] = (('event_index'), PeakStart[:, i])
    my_ds['PkStart_ch' + str(i)].attrs["long_name"] = "Peak start for channel %d" % i
    my_ds['PkStart_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['GaussChiSq_ch' + str(i)] = (('event_index'), GaussChiSq[:, i])
    my_ds['GaussChiSq_ch' + str(i)].attrs["long_name"] = "Chisquare value for channel %d" % i
    my_ds['GaussChiSq_ch' + str(i)].attrs["_FillValue"] = np.nan

    # Incandescence stuff
    where_valid = np.logical_and(my_ds['PkPos_ch2'] != np.nan,
                                 my_ds['PkPos_ch1'] != np.nan)
    my_ds['IncanPkOffsetch1ch2'] = my_ds['PkPos_ch2'] - my_ds['PkPos_ch1']
    my_ds['IncanPkOffsetch1ch2'][~where_valid] = np.nan
    where_valid = np.logical_and(my_ds['PkPos_ch5'] != np.nan,
                                 my_ds['PkPos_ch6'] != np.nan)
    my_ds['IncanPkOffsetch5ch6'] = my_ds['PkPos_ch6'] - my_ds['PkPos_ch5']
    my_ds['IncanPkOffsetch5ch6'][~where_valid] = np.nan
    my_ds['IncanRatioch1ch2'] = (('event_index'), IncanRatioch1ch2)
    my_ds['IncanRatioch1ch2'].attrs["long_name"] = "Incandescence ratio ch1, ch2"
    my_ds['IncanRatioch1ch2'].attrs["_FillValue"] = np.nan
    my_ds['IncanRatioch5ch6'] = (('event_index'), IncanRatioch5ch6)
    my_ds['IncanRatioch5ch6'].attrs["long_name"] = "Incandescence ratio ch5, ch6"
    my_ds['IncanRatioch5ch6'].attrs["_FillValue"] = np.nan

    # First do initial filter step
    scat_reject = np.logical_and.reduce(
        (np.isfinite(my_ds['PkHt_ch0'].values), np.isfinite(my_ds['PkFWHM_ch0'].values),
         np.isfinite(my_ds['PkPos_ch0'].values)))
    scat_reject = np.logical_and.reduce(
         (scat_reject, np.isfinite(my_ds['PkHt_ch3'].values),
         np.isfinite(my_ds['PkFWHM_ch3'].values), np.isfinite(my_ds['PkSplitPos_ch3'].values)))
    incan_reject = np.logical_and.reduce(
        (np.isfinite(my_ds['PkHt_ch1'].values), np.isfinite(my_ds['PkEnd_ch1'].values),
         np.isfinite(my_ds['PkStart_ch1'].values), np.isfinite(my_ds['PkPos_ch1'].values)))
    incan_reject = np.logical_and.reduce(
        (incan_reject, np.isfinite(my_ds['PkHt_ch2'].values), np.isfinite(my_ds['PkEnd_ch2'].values),
         np.isfinite(my_ds['PkStart_ch2'].values), np.isfinite(my_ds['PkPos_ch2'].values)))

    scat_reject_key = np.where(~scat_reject, 1, 0)
    incan_reject_key = np.where(~incan_reject, 1, 0)
    DMTglobals = DMTGlobals()
    # Then we apply criteria to max.min peak heights
    scat_reject_reason2 = np.logical_and.reduce((~scat_reject, my_ds['PkHt_ch0'].values < DMTglobals.ScatMinPeakHt,
                                                 my_ds['PkHt_ch3'].values < DMTglobals.ScatMinPeakHt))
    scat_reject_reason3 = np.logical_and.reduce((~scat_reject_reason2, my_ds['PkHt_ch0'].values > DMTglobals.ScatMaxPeakHt,
                                                 my_ds['PkHt_ch3'].values > DMTglobals.ScatMaxPeakHt))
    scat_reject_reason4 = np.logical_and.reduce((~scat_reject_reason3, my_ds['PkFWHM_ch0'].values < DMTglobals.ScatMinWidth,
                                                 my_ds['PkFWHM_ch3'].values < DMTglobals.ScatMinWidth))
    scat_reject_reason5 = np.logical_and.reduce((~scat_reject_reason4, my_ds['PkFWHM_ch0'].values < DMTglobals.ScatMaxWidth,
                                                 my_ds['PkFWHM_ch3'].values > DMTglobals.ScatMaxWidth))
    scat_reject_reason6 = np.logical_and.reduce((~scat_reject_reason5, my_ds['PkPos_ch0'].values < DMTglobals.ScatMinPeakPos,
                                                 my_ds['PkPos_ch3'].values < DMTglobals.ScatMinPeakPos))
    scat_reject_reason7 = np.logical_and.reduce((~scat_reject_reason6, my_ds['PkPos_ch0'].values < DMTglobals.ScatMaxPeakPos,
                                                 my_ds['PkPos_ch3'].values > DMTglobals.ScatMaxPeakPos))

    incan_reject_reason2 = np.logical_and.reduce((~incan_reject, my_ds['PkHt_ch1'].values < DMTglobals.IncanMinPeakHt,
                                                 my_ds['PkHt_ch2'].values < DMTglobals.IncanMinPeakHt))
    incan_reject_reason3 = np.logical_and.reduce((~incan_reject_reason2, my_ds['PkHt_ch1'].values > DMTglobals.IncanMaxPeakHt,
                                                 my_ds['PkHt_ch2'].values > DMTglobals.IncanMaxPeakHt))
    width1 = my_ds['PkEnd_ch1'].values-my_ds['PkStart_ch1'].values
    width2 = my_ds['PkEnd_ch2'].values-my_ds['PkStart_ch2'].values
    incan_reject_reason4 = np.logical_and.reduce((~incan_reject_reason3, width1 < DMTglobals.IncanMinWidth,
                                                 width2 < DMTglobals.IncanMinWidth))
    incan_reject_reason5 = np.logical_and.reduce((~incan_reject_reason4, width1 < DMTglobals.IncanMaxWidth,
                                                 width2  > DMTglobals.IncanMaxWidth))
    incan_reject_reason6 = np.logical_and.reduce(
        (~incan_reject_reason5, my_ds['PkPos_ch1'].values < DMTglobals.IncanMinPeakPos,
         my_ds['PkPos_ch2'].values < DMTglobals.IncanMinPeakPos))
    incan_reject_reason7 = np.logical_and.reduce(
        (~incan_reject_reason6, my_ds['PkPos_ch1'].values < DMTglobals.IncanMaxPeakPos,
         my_ds['PkPos_ch2'].values > DMTglobals.IncanMaxPeakPos))
    incan_reject_reason8 = np.logical_and.reduce(
        (~incan_reject_reason7, my_ds['IncanRatioch1ch2'].values < DMTglobals.IncanMinPeakRatio,
         my_ds['IncanRatioch5ch6'].values > DMTglobals.IncanMinPeakRatio))
    incan_reject_reason9 = np.logical_and.reduce(
        (~incan_reject_reason8, my_ds['IncanRatioch1ch2'].values > DMTglobals.IncanMaxPeakRatio,
         my_ds['IncanRatioch5ch6'].values > DMTglobals.IncanMaxPeakRatio))
    incan_reject_reason10 = np.logical_and.reduce(
        (~incan_reject_reason9, my_ds['IncanPkOffsetch1ch2'].values > DMTglobals.IncanMaxPeakOffset,
         my_ds['IncanPkOffsetch5ch6'].values > DMTglobals.IncanMaxPeakOffset))

    scat_reject_key[scat_reject_reason2] = 2
    scat_reject_key[scat_reject_reason3] = 3
    scat_reject_key[scat_reject_reason4] = 4
    scat_reject_key[scat_reject_reason5] = 5
    scat_reject_key[scat_reject_reason6] = 6
    scat_reject_key[scat_reject_reason7] = 7
    incan_reject_key[incan_reject_reason2] = 2
    incan_reject_key[incan_reject_reason3] = 3
    incan_reject_key[incan_reject_reason4] = 4
    incan_reject_key[incan_reject_reason5] = 5
    incan_reject_key[incan_reject_reason6] = 6
    incan_reject_key[incan_reject_reason7] = 7
    incan_reject_key[incan_reject_reason8] = 8
    incan_reject_key[incan_reject_reason9] = 9
    incan_reject_key[incan_reject_reason10] = 10

    my_ds['ScatRejectKey'] = (('event_index'), scat_reject_key)
    my_ds['ScatRejectKey'].attrs["long_name"] = "Scattering reject flag"
    my_ds['ScatRejectKey'].attrs["_FillValue"] = np.nan
    my_ds['IncanRejectKey'] = (('event_index'), incan_reject_key)
    my_ds['IncanRejectKey'].attrs["long_name"] = "Incandescence reject flag"
    my_ds['IncanRejectKey'].attrs["_FillValue"] = np.nan

    print(str(num_records) + ' records processed in ' + str(time.time()-start_time) + ' s')
    return my_ds


def _gaus(x, a, x0, sigma, base):
    return a * np.exp(-((x - x0)/sigma)**2) + base


def _fit_record_gaussian(my_ds, record_number):
    """ Only used for channel 0."""
    bins = np.arange(0, 100, 1.)
    p0 = []
    chn = 0
    data = my_ds['Data_ch' + str(chn)].isel(event_index=record_number).values
    height = np.nan
    pos = np.nan
    start = np.nan
    error = np.nan
    try:
        valid_inds = np.where(data < 32767)
        data_fit = data[valid_inds]
        bins_fit = bins[valid_inds]
        p0 = [data_fit.max()-data_fit.min(), float(np.argmax(data_fit)), 20., np.nanmean(data_fit)]
        coeff, var_matrix = curve_fit(_gaus, bins_fit, data_fit, p0=p0)
        amplitude = coeff[0]
        peakpos = coeff[1]
        width = coeff[2]*(2.35482/np.sqrt(2))
        base = coeff[3]
        fit_data = _gaus(np.array(bins, dtype=np.float64), *coeff)
        chi2 = chisquare(np.array(data, dtype='float64'), f_exp=np.array(fit_data, dtype='float64'))
        if not (amplitude > 1 and peakpos > 0 and
                peakpos < len(data) and width < len(data) and width < amplitude):
            amplitude = np.nan
            width = np.nan
            peakpos = np.nan
            base = np.nan
            chi2 = np.nan
    except RuntimeError as e:
        amplitude = np.nan
        width = np.nan
        peakpos = np.nan
        base = np.nan
        chi2 = np.nan
        error = 1
    except MemoryError:
        print(e)
        amplitude = np.nan
        width = np.nan
        peakpos = np.nan
        base = np.nan
        chi2 = np.nan
        error = 2

    if np.isfinite(base):
        try:
            height = data.max() - base
            pos = np.argmax(data)
        except ValueError:
            height = np.nan
            pos = np.nan

        try:
            start = np.where((data - base) < 50)[0]
            if start == []:
                start = np.nan
            else:
                start = start[start <= pos][-1]
        except IndexError:
            start = np.nan


        # Filter out bad points
        bad = ~np.logical_and.reduce(
            (height > 1, peakpos > 0, peakpos < len(bins)))
        if bad:
            amplitude = np.nan
            base = np.nan
            peakpos = np.nan
            width = np.nan
            chi2 = np.nan
            pos = np.nan
    else:
        height = np.nan
        pos = np.nan
        start = np.nan

    fit_coeffs = {'amplitude': amplitude, 'peakpos': peakpos,
                  'width': width, 'base': base, 'chi2': chi2,
                  'height': height, 'pos': pos, 'start': start,
                  'error': error}
    return fit_coeffs


def _fit_record_incan_ave_base(my_ds, record_number, channel, num_trig_pts):
    """ Channels 1, 2, 6, 7"""
    num_base_pts_2_avg_backup = 20
    num_pts = len(my_ds['Data_ch' + str(channel)].values)
    if num_trig_pts != -1:
        num_base_pts_2_avg = round(0.8*num_trig_pts)
        if((not np.isfinite(num_base_pts_2_avg)) or
                num_base_pts_2_avg <= 0 or
            num_base_pts_2_avg >= num_pts):
            num_base_pts_2_avg = num_base_pts_2_avg_backup
    else:
        num_base_pts_2_avg = num_base_pts_2_avg_backup

    base = np.nan
    height = np.nan
    pos = np.nan
    start = np.nan
    end = np.nan
    half_rise = np.nan
    half_decay = np.nan
    peak2area = np.nan
    data = my_ds['Data_ch' + str(channel)].isel(event_index=record_number).values
    base = np.nanmean(data[0:num_base_pts_2_avg])
    data2 = data + abs(base)
    V_max = data.max()
    V_maxloc = np.where(data == V_max)[0][0]
    try:
        peak2area = np.max(data2)/np.sum(data2[20:81])
    except ZeroDivisionError:
        peak2area = np.nan

    if channel in [1, 5]:
        upperbound = 100000.
        lowerbound = 0.063
    else:
        upperbound = 1000000.
        lowerbound = -1000000.

    #if(len(V_maxloc) > 1):

    if((V_max - base) > 1 and V_maxloc > 0 and V_maxloc < len(data) and
       peak2area > lowerbound and peak2area < upperbound):
        height = V_max - base
        pos = V_maxloc
    else:
        base = np.nan
        height = np.nan
        pos = np.nan
        start = np.nan
        end = np.nan
        half_rise = np.nan
        half_decay = np.nan
        peak2area = np.nan

    if np.isfinite(pos):
        diffs = data - base
        under_peak = np.where(diffs < 5)[0]
        under_half = np.where(diffs <= height*0.5)[0]
        try:
            start = under_peak[under_peak <= pos][-1]
        except IndexError:
            start = np.nan
        try:
            end = under_peak[under_peak >= pos][0]
        except IndexError:
            end = np.nan
        try:
            half_rise = under_half[under_half <= pos][-1]
        except IndexError:
            half_rise = np.nan

        try:
            half_decay = under_half[under_half >= pos][0]
        except IndexError:
            half_decay = np.nan

    fit_coeffs = {'base': base, 'height': height, 'pos': pos, 'start': start,
                  'end': end, 'half_rise': half_rise, 'half_decay': half_decay,
                  'peak2area': peak2area}

    return fit_coeffs


def _split_scatter_fit(my_ds, record_number, channel):
    """ Used for channels 3, 7"""
    base = np.nan
    height = np.nan
    peak = np.nan
    start = np.nan
    pos = np.nan
    num_base_pts_2_avg = 20
    data = my_ds['Data_ch' + str(channel)].isel(event_index=record_number).values
    V_maxloc = np.argmax(data)
    V_minloc = np.argmin(data)
    if V_maxloc < V_minloc:
        data = -data
    base = np.nanmean(data[0:num_base_pts_2_avg])
    V_max = data.max()
    if((V_max - base) > 1 and V_maxloc < len(data)):
        height = V_max - base
        pos = np.argmax(data)
    else:
        base = np.nan
        pos = np.nan
        start = np.nan
        height = np.nan

    data = data - base
    data = np.abs(data) + data
    if np.isfinite(pos):
        lt5 = np.where(data < 5)[0]
        try:
            start = lt5[lt5 <= pos][-1]
        except IndexError:
            start = np.nan
    fit_coeffs = {'base': base, 'height': height, 'pos': pos, 'start': start}

    return fit_coeffs

def _gaussian_sat_fit(my_ds, record_number):
    channel = 4
    base = np.nan
    fitamplitude = np.nan
    fitpos = np.nan
    height = np.nan
    pos = np.nan
    chi2 = np.nan
    width = np.nan
    error = np.nan
    error_thrown = False
    start = np.nan
    clipped_wave = False
    global_vars = DMTGlobals()
    data = my_ds['Data_ch' + str(channel)].isel(event_index=record_number).values

    if data.max() - data.min() >= global_vars.ScatMaxPeakHt:
        temp1 = data.astype(float)
        temp1[temp1 == data.max()] = np.nan
        clipped_wave = True

    bins = np.arange(0, 100, 1.)
    try:
        p0 = [data.max()-data.min(), np.argmax(data), np.argmax(data), np.nanmin(data)]
        coeff, var_matrix = curve_fit(_gaus, bins, data, p0=p0)
        if clipped_wave:
            p0[1] = coeff[1]
            bins_fit = bins[np.isfinite(temp1)]
            temp1_fit = temp1[np.isfinite(temp1)]

            bounds = ([0, 0, 0, -np.inf], [np.inf, len(data), len(data), np.inf])
            coeff, var_matrix = curve_fit(_gaus, bins_fit, temp1_fit, p0=p0)

        fit_data = _gaus(bins, *coeff)
        chi2 = chisquare(np.array(data, dtype='float64'), f_exp=np.array(fit_data, dtype='float64'))
        fitamplitude = coeff[0]
        fitpos = coeff[1]
        width = coeff[2]*(2.35482/np.sqrt(2))
        base = coeff[3]
        if not (fitamplitude > 1 and fitpos > 0 and width < len(data) and width < fitamplitude):
            fitamplitude = np.nan
            fitpos = np.nan
            width = np.nan
            chi2 = np.nan
            base = np.nan
    except RuntimeError as e:
        error = 1
        error_thrown = True
    except MemoryError:
        error = 2

    if np.isfinite(base):
        height = data.max() - base
        pos = np.argmax(data)
        if not (height > 1 and pos > 0 and pos < len(data)):
            height = np.nan
            pos = np.nan

        lt5 = np.where(data - base < 50)[0]
        try:
            start = lt5[lt5 <= pos][-1]
        except IndexError:
            start = np.nan

    fit_coeffs = {'base': base, 'fitamplitude': fitamplitude, 'fitpos': fitpos, 'start': start,
                  'pos': pos, 'chi2': chi2, 'error_thrown': error_thrown, 'width': width,
                  'height': height, 'error': error}

    return fit_coeffs





