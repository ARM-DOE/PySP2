import numpy as np
import time
import dask.bag as db

from scipy.optimize import curve_fit
from .DMTGlobals import DMTGlobals


def _do_fit_records(my_ds, i, num_trig_pts, debug=True):
    if debug and i % 1000 == 0:
        print("Processing record %d" % i)

    FtAmp = np.zeros(2)
    FtPos = np.zeros(2)
    Base = np.zeros(2)
    PeakHeight = np.zeros(2)
    PeakPos = np.zeros(2)
    GaussChiSq = np.zeros(2)
    PeakStart = np.zeros(2)
    PeakEnd = np.zeros(2)
    Width = np.zeros(2)
    HalfRise = np.zeros(2)
    HalfDecay = np.zeros(2)
    Peak2Area = np.zeros(2)
    Error = np.zeros(2)

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

    # Channel 4 is a special snowflake
    coeffs = _gaussian_sat_fit(my_ds, i)
    FtAmp[1] = coeffs['fitamplitude']
    Base[1] = coeffs['base']
    PeakHeight[1] = coeffs['height']
    PeakPos[1] = coeffs['pos']
    FtPos[1] = coeffs['fitpos']
    PeakStart[1] = coeffs['start']
    GaussChiSq[1] = coeffs['chi2']
    Width[1] = coeffs['width']
    Error[1] = coeffs['error']

    return (FtAmp, FtPos, Base, PeakHeight, PeakPos, GaussChiSq, PeakStart,
            PeakEnd, Width, HalfRise, HalfDecay, Peak2Area, Error)


def _calc_incan_ratio(my_ds, ch1, ch2):
    data_ch1 = my_ds['Data_ch' + str(ch1)].values
    data_ch2 = my_ds['Data_ch' + str(ch2)].values
    PeakPos_ch1 = my_ds['PkPos_ch%d' % ch1].values
    halfDecay_ch1 = my_ds['PkHalfDecay_ch%d' % ch1].values
    PeakPos_ch1_tile = np.tile(PeakPos_ch1, (data_ch1.shape[1], 1)).T
    halfDecay_ch1_tile = np.tile(halfDecay_ch1, (data_ch1.shape[1], 1)).T
    Base_ch1 = my_ds['Base_ch%d' % ch1].values
    Base_ch2 = my_ds['Base_ch%d' % ch2].values
    Base_ch1_tile = np.tile(Base_ch1, (data_ch1.shape[1], 1)).T
    Base_ch2_tile = np.tile(Base_ch2, (data_ch1.shape[1], 1)).T
    finite_mask = np.logical_and(
        np.isfinite(PeakPos_ch1_tile), halfDecay_ch1_tile)
    finite_mask = np.logical_and(finite_mask, data_ch2 - Base_ch2_tile > 0)
    counting_up = np.tile(np.arange(data_ch1.shape[1]), (data_ch1.shape[0], 1))
    range_mask = np.logical_and(
        counting_up >= PeakPos_ch1_tile, counting_up <= halfDecay_ch1_tile)
    data_ch2 = np.where(
        np.logical_and(finite_mask, range_mask), data_ch2, np.nan)
    data_ch1 = np.where(
        np.logical_and(finite_mask, range_mask), data_ch1, np.nan)
    ratio = np.nanmean(
        (data_ch1 - Base_ch1_tile)/(data_ch2 - Base_ch2_tile), axis=1)
    return ratio


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

    for i in [3, 7]:
        coeffs = _split_scatter_fit(my_ds, i)
        Base2 = coeffs['base']
        PeakHeight2 = coeffs['height']
        PeakPos2 = coeffs['pos']
        PeakStart2 = coeffs['start']
        my_ds['Base_ch' + str(i)] = (('event_index'), Base2)
        my_ds['Base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
        my_ds['Base_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkHt_ch' + str(i)] = (('event_index'), PeakHeight2)
        my_ds['PkHt_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
        my_ds['PkHt_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkSplitPos_ch' + str(i)] = (('event_index'), PeakStart2)
        my_ds['PkSplitPos_ch' + str(i)].attrs["long_name"] = "Peak start position for channel %d" % i
        my_ds['PkSplitPos_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkPos_ch' + str(i)] = (('event_index'), PeakPos2)
        my_ds['PkPos_ch' + str(i)].attrs["long_name"] = "Peak split position for channel %d" % i
        my_ds['PkPos_ch' + str(i)].attrs["_FillValue"] = np.nan

    for i in [1, 2, 5, 6]:
        coeffs = _fit_record_incan_ave_base(my_ds, i, num_trig_pts)
        Base = coeffs['base']
        PeakHeight2 = coeffs['height']
        PeakPos2 = coeffs['pos']
        PeakStart2 = coeffs['start']
        PeakEnd2 = coeffs['end']
        HalfRise2 = coeffs['half_rise']
        HalfDecay2 = coeffs['half_decay']
        Peak2Area2 = coeffs['peak2area']

        my_ds['Base_ch' + str(i)] = (('event_index'), Base)
        my_ds['Base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
        my_ds['Base_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkHt_ch' + str(i)] = (('event_index'), PeakHeight2)
        my_ds['PkHt_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
        my_ds['PkHt_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkHalfRise_ch' + str(i)] = (('event_index'), HalfRise2)
        my_ds['PkHalfRise_ch' + str(i)].attrs["long_name"] = "Point where rise is at 1/2 height for channel %d" % i
        my_ds['PkHalfRise_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['Peak2area_ch' + str(i)] = (('event_index'), Peak2Area2)
        my_ds['Peak2area_ch' + str(i)].attrs["long_name"] = "Peak 2 area for channel %d" % i
        my_ds['Peak2area_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkHalfDecay_ch' + str(i)] = (('event_index'), HalfDecay2)
        my_ds['PkHalfDecay_ch' + str(i)].attrs["long_name"] = "Point where decay is at 1/2 height for channel %d" % i
        my_ds['PkHalfDecay_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkPos_ch' + str(i)] = (('event_index'), PeakPos2)
        my_ds['PkPos_ch' + str(i)].attrs["long_name"] = "Peak position for channel %d" % i
        my_ds['PkPos_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkStart_ch' + str(i)] = (('event_index'), PeakStart2)
        my_ds['PkStart_ch' + str(i)].attrs["long_name"] = "Peak start for channel %d" % i
        my_ds['PkStart_ch' + str(i)].attrs["_FillValue"] = np.nan
        my_ds['PkEnd_ch' + str(i)] = (('event_index'), PeakEnd2)
        my_ds['PkEnd_ch' + str(i)].attrs["long_name"] = "Peak end for channel %d" % i
        my_ds['PkEnd_ch' + str(i)].attrs["_FillValue"] = np.nan

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
    Width = np.stack([x[8] for x in proc_records])

    # Channel 0
    i = 0
    my_ds['FtAmp_ch' + str(i)] = (('event_index'), FtAmp[:, 0])
    my_ds['FtAmp_ch' + str(i)].attrs["long_name"] = "Fit Amplitude for channel %d" % i
    my_ds['FtAmp_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['FtPos_ch' + str(i)] = (('event_index'), FtPos[:, 0])
    my_ds['FtPos_ch' + str(i)].attrs["long_name"] = "Fit Position for channel %d" % i
    my_ds['FtPos_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['Base_ch' + str(i)] = (('event_index'), Base[:, 0])
    my_ds['Base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
    my_ds['Base_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkHt_ch' + str(i)] = (('event_index'), PeakHeight[:, 0])
    my_ds['PkHt_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
    my_ds['PkHt_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkFWHM_ch' + str(i)] = (('event_index'), Width[:, 0])
    my_ds['PkFWHM_ch' + str(i)].attrs["long_name"] = "Width for channel %d" % i
    my_ds['PkFWHM_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkPos_ch' + str(i)] = (('event_index'), PeakPos[:, 0])
    my_ds['PkPos_ch' + str(i)].attrs["long_name"] = "Peak position for channel %d" % i
    my_ds['PkPos_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkStart_ch' + str(i)] = (('event_index'), PeakStart[:, 0])
    my_ds['PkStart_ch' + str(i)].attrs["long_name"] = "Peak start for channel %d" % i
    my_ds['PkStart_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['GaussChiSq_ch' + str(i)] = (('event_index'), GaussChiSq[:, 0])
    my_ds['GaussChiSq_ch' + str(i)].attrs["long_name"] = "Chisquare value for channel %d" % i
    my_ds['GaussChiSq_ch' + str(i)].attrs["_FillValue"] = np.nan

    # Channel 4
    i = 4
    my_ds['FtAmp_ch' + str(i)] = (('event_index'), FtAmp[:, 1])
    my_ds['FtAmp_ch' + str(i)].attrs["long_name"] = "Amplitude for channel %d" % i
    my_ds['FtAmp_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['Base_ch' + str(i)] = (('event_index'), Base[:, 1])
    my_ds['Base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
    my_ds['Base_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkHt_ch' + str(i)] = (('event_index'), PeakHeight[:, 1])
    my_ds['PkHt_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
    my_ds['PkHt_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkFWHM_ch' + str(i)] = (('event_index'), Width[:, 1])
    my_ds['PkFWHM_ch' + str(i)].attrs["long_name"] = "Width for channel %d" % i
    my_ds['PkFWHM_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkPos_ch' + str(i)] = (('event_index'), PeakPos[:, 1])
    my_ds['PkPos_ch' + str(i)].attrs["long_name"] = "Peak position for channel %d" % i
    my_ds['PkPos_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['FtPos_ch' + str(i)] = (('event_index'), FtPos[:, 1])
    my_ds['FtPos_ch' + str(i)].attrs["long_name"] = "Fit position for channel %d" % i
    my_ds['FtPos_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['PkStart_ch' + str(i)] = (('event_index'), PeakStart[:, 1])
    my_ds['PkStart_ch' + str(i)].attrs["long_name"] = "Peak start for channel %d" % i
    my_ds['PkStart_ch' + str(i)].attrs["_FillValue"] = np.nan
    my_ds['GaussChiSq_ch' + str(i)] = (('event_index'), GaussChiSq[:, 1])
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

    IncanRatioch1ch2 = _calc_incan_ratio(my_ds, 1, 2)
    my_ds['IncanRatioch1ch2'] = (('event_index'), IncanRatioch1ch2)
    my_ds['IncanRatioch1ch2'].attrs["long_name"] = "Incandescence ratio ch1, ch2"
    my_ds['IncanRatioch1ch2'].attrs["_FillValue"] = np.nan
    IncanRatioch5ch6 = _calc_incan_ratio(my_ds, 5, 6)
    my_ds['IncanRatioch5ch6'] = (('event_index'), IncanRatioch5ch6)
    my_ds['IncanRatioch5ch6'].attrs["long_name"] = "Incandescence ratio ch5, ch6"
    my_ds['IncanRatioch5ch6'].attrs["_FillValue"] = np.nan

    # First do initial filter step
    scat_reject = np.logical_or.reduce(
        (~np.isfinite(my_ds['PkHt_ch0'].values), ~np.isfinite(my_ds['PkFWHM_ch0'].values),
         ~np.isfinite(my_ds['PkPos_ch0'].values)))

    incan_reject = np.logical_or.reduce(
        (~np.isfinite(my_ds['PkHt_ch1'].values), ~np.isfinite(my_ds['PkEnd_ch1'].values),
         ~np.isfinite(my_ds['PkStart_ch1'].values), ~np.isfinite(my_ds['PkPos_ch1'].values, ~np.isfinite(my_ds['IncanRatioch1ch2'].values))))


    scat_reject_key = np.where(scat_reject, 1, 0)
    incan_reject_key = np.where(incan_reject, 1, 0)
    DMTglobals = DMTGlobals()
    # Then we apply criteria to max.min peak heights
    scat_reject_reason2 = np.logical_and.reduce((~scat_reject, my_ds['PkHt_ch0'].values < DMTglobals.ScatMinPeakHt1))
    prev_scat_reject = np.logical_or(scat_reject, scat_reject_reason2)
    scat_reject_reason3 = np.logical_and.reduce((~prev_scat_reject, my_ds['PkHt_ch0'].values > DMTglobals.ScatMaxPeakHt1))
    prev_scat_reject = np.logical_or(prev_scat_reject, scat_reject_reason3)
    scat_reject_reason4 = np.logical_and(~prev_scat_reject, my_ds['PkFWHM_ch0'].values < DMTglobals.ScatMinWidth)
    prev_scat_reject = np.logical_or(prev_scat_reject, scat_reject_reason4)
    scat_reject_reason5 = np.logical_and(~prev_scat_reject, my_ds['PkFWHM_ch0'].values > DMTglobals.ScatMaxWidth)
    prev_scat_reject = np.logical_or(prev_scat_reject, scat_reject_reason5)
    scat_reject_reason6 = np.logical_and.reduce((~prev_scat_reject,
                                                 my_ds['PkPos_ch0'].values < DMTglobals.ScatMinPeakPos))
    prev_scat_reject = np.logical_or(prev_scat_reject, scat_reject_reason6)
    scat_reject_reason7 = np.logical_and.reduce((~prev_scat_reject,
                                                 my_ds['PkPos_ch0'].values > DMTglobals.ScatMaxPeakPos))

    incan_reject_reason2 = np.logical_and(
        ~incan_reject, my_ds['PkHt_ch1'].values < DMTglobals.IncanMinPeakHt1)
    prev_incan_reject = np.logical_or(incan_reject, incan_reject_reason2)
    incan_reject_reason3 = np.logical_and(
        ~prev_incan_reject, my_ds['PkHt_ch1'].values > DMTglobals.IncanMaxPeakHt1)
    width1 = my_ds['PkEnd_ch1'].values - my_ds['PkStart_ch1'].values

    prev_incan_reject = np.logical_or(prev_incan_reject, incan_reject_reason3)
    incan_reject_reason4 = np.logical_and.reduce((~prev_incan_reject, width1 < DMTglobals.IncanMinWidth))
    prev_incan_reject = np.logical_or(prev_incan_reject, incan_reject_reason4)
    incan_reject_reason5 = np.logical_and.reduce((~prev_incan_reject, width1 > DMTglobals.IncanMaxWidth))
    prev_incan_reject = np.logical_or(prev_incan_reject, incan_reject_reason5)
    incan_reject_reason6 = np.logical_and.reduce(
        (~prev_incan_reject, my_ds['PkPos_ch1'].values < DMTglobals.IncanMinPeakPos,
         ))
    prev_incan_reject = np.logical_or(prev_incan_reject, incan_reject_reason6)
    incan_reject_reason7 = np.logical_and.reduce(
        (~prev_incan_reject, my_ds['PkPos_ch1'].values > DMTglobals.IncanMaxPeakPos,
         ))
    prev_incan_reject = np.logical_or(prev_incan_reject, incan_reject_reason7)
    incan_reject_reason8 = np.logical_and.reduce(
        (~prev_incan_reject, np.logical_or(my_ds['IncanRatioch1ch2'].values < DMTglobals.IncanMinPeakRatio,
         my_ds['IncanRatioch5ch6'].values < DMTglobals.IncanMinPeakRatio)))
    prev_incan_reject = np.logical_or(prev_incan_reject, incan_reject_reason8)
    incan_reject_reason9 = np.logical_and.reduce(
        (~prev_incan_reject, np.logical_or(my_ds['IncanRatioch1ch2'].values > DMTglobals.IncanMaxPeakRatio,
         my_ds['IncanRatioch5ch6'].values > DMTglobals.IncanMaxPeakRatio)))
    prev_incan_reject = np.logical_or(prev_incan_reject, incan_reject_reason9)
    incan_reject_reason10 = np.logical_and.reduce(
        (~prev_incan_reject, np.logical_or(my_ds['IncanPkOffsetch1ch2'].values > DMTglobals.IncanMaxPeakOffset,
         my_ds['IncanPkOffsetch5ch6'].values > DMTglobals.IncanMaxPeakOffset)))
    prev_incan_reject = np.logical_or(prev_incan_reject, incan_reject_reason10)

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
    return a * np.exp(-((x - x0)**2/(2 * sigma**2))) + base


def _fit_record_gaussian(my_ds, record_number):
    """ Only used for channel 0."""
    bins = np.arange(0, 100, 1.)
    p0 = []
    chn = 0
    data = my_ds['Data_ch' + str(chn)].values[record_number]
    height = np.nan
    pos = np.nan
    start = np.nan
    error = np.nan
    try:
        data_fit = data
        bins_fit = bins
        p0 = np.array([data_fit.max()-data_fit.min(), np.argmax(data_fit), 20., np.nanmin(data_fit)]).astype(float)

        coeff, var_matrix = curve_fit(_gaus, bins_fit, data_fit, p0=p0, method='lm', maxfev=40, ftol=1e-3)
        amplitude = coeff[0]
        peakpos = coeff[1]
        width = coeff[2]*(2.35482)
        base = coeff[3]
        fit_data = _gaus(np.array(bins, dtype=np.float64), *coeff)
        chi2 = chisquare(np.array(data, dtype='float64'), f_exp=np.array(fit_data, dtype='float64'))
        if not (amplitude > 1 and peakpos > 0 and
                peakpos < len(data) and width < len(data) and width < amplitude and width > 0):
            amplitude = np.nan
            width = np.nan
            peakpos = np.nan
            base = np.nan
            chi2 = np.nan
    except RuntimeError:
        amplitude = np.nan
        width = np.nan
        peakpos = np.nan
        base = np.nan
        chi2 = np.nan
        error = 1
    except MemoryError:
        amplitude = np.nan
        width = np.nan
        peakpos = np.nan
        base = np.nan
        chi2 = np.nan
        error = 2

    if np.isfinite(base):
        try:
            height = data_fit.max() - base
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


def _fit_record_incan_ave_base(my_ds, channel, num_trig_pts):
    """ Channels 1, 2, 6, 7"""
    num_base_pts_2_avg_backup = 20
    num_pts = my_ds['Data_ch' + str(channel)].values.shape[1]
    if num_trig_pts != -1:
        num_base_pts_2_avg = round(0.8*num_trig_pts)
        if((not np.isfinite(num_base_pts_2_avg)) or
            num_base_pts_2_avg <= 0 or num_base_pts_2_avg >= num_pts):
            num_base_pts_2_avg = num_base_pts_2_avg_backup
    else:
        num_base_pts_2_avg = num_base_pts_2_avg_backup

    data = my_ds['Data_ch' + str(channel)].values.astype(int)
    base = np.mean(data[:, 0:num_base_pts_2_avg], axis=1)
    data2 = data + abs(np.tile(base, (data.shape[1], 1))).T
    V_max = data.max(axis=1)
    V_maxloc = np.argmax(data, axis=1)
    denominator = np.sum(data2[:, 20:81], axis=1)
    peak2area = np.max(data2, axis=1)/denominator

    conditions = np.logical_and.reduce(
        (V_max - base > 1, V_maxloc > 0, V_maxloc < data.shape[1]))
    height = np.where(conditions, V_max - base, np.nan)
    pos = np.where(conditions, V_maxloc, np.nan)
    base = np.where(conditions, base, np.nan)

    diffs = data - np.tile(base, (data.shape[1], 1)).T
    pos_tile = np.tile(pos, (data.shape[1], 1)).T
    height_tile = np.tile(height, (data.shape[1], 1)).T
    counting_up = np.tile(np.arange(data.shape[1]), (data.shape[0], 1))
    start_numbers = np.where(np.logical_and(diffs < 5, counting_up <= pos_tile), counting_up, -1)
    start = start_numbers.max(axis=1).astype(float)
    start[start == -1] = np.nan
    end_numbers = np.where(np.logical_and(diffs < 5, counting_up >= pos_tile), counting_up, 9999)
    end = end_numbers.min(axis=1).astype(float)
    end[end == 9999] = np.nan
    start_numbers = np.where(np.logical_and(diffs <= 0.5*height_tile, counting_up <= pos_tile), counting_up, -1)
    half_rise = start_numbers.max(axis=1).astype(float)
    half_rise[half_rise == -1] = np.nan
    end_numbers = np.where(np.logical_and(diffs <= 0.5*height_tile, counting_up >= pos_tile), counting_up, 9999)
    half_decay = end_numbers.min(axis=1).astype(float)
    half_decay[half_decay == 9999] = np.nan
    start = np.where(conditions, start, np.nan)
    end = np.where(conditions, end, np.nan)
    half_rise = np.where(conditions, half_rise, np.nan)
    half_decay = np.where(conditions, half_decay, np.nan)

    fit_coeffs = {'base': base, 'height': height, 'pos': pos, 'start': start,
                  'end': end, 'half_rise': half_rise, 'half_decay': half_decay,
                  'peak2area': peak2area}

    return fit_coeffs


def _split_scatter_fit(my_ds, channel):
    """ Used for channels 3, 7"""

    num_base_pts_2_avg = 20
    data = my_ds['Data_ch' + str(channel)].values
    V_maxloc = np.argmax(data, axis=1)
    V_minloc = np.argmin(data, axis=1)
    data[V_maxloc < V_minloc, :] = -data[V_maxloc < V_minloc, :]
    base = np.nanmean(data[:, 0:num_base_pts_2_avg], axis=1)
    V_max = data.max(axis=1)
    conditions = np.logical_and.reduce(((V_max - base) > 1, V_maxloc < len(data), V_maxloc > 0))
    height = np.where(conditions, V_max - base, np.nan)
    pos = np.where(conditions, np.argmax(data, axis=1), np.nan)
    start = np.zeros_like(height)
    start[~conditions] = np.nan
    height[~conditions] = np.nan
    data = data - np.tile(base, (data.shape[1], 1)).T
    counting_up = np.tile(np.arange(data.shape[1]), (data.shape[0], 1))
    data = np.abs(data) + data
    pos_tile = np.tile(pos, (data.shape[1], 1)).T
    counting_up = np.where(np.logical_and(data < 5, counting_up <= pos_tile), counting_up, -1)
    start = counting_up.max(axis=1)
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
    data = my_ds['Data_ch' + str(channel)].values[record_number]

    if data.max() - data.min() >= global_vars.ScatMaxPeakHt1:
        temp1 = data.astype(float)
        temp1[temp1 == data.max()] = np.nan
        clipped_wave = True
    else:
        temp1 = data.astype(float)
        clipped_wave = False

    bins = np.arange(0, 100, 1.)
    try:
        bins_fit = bins[np.isfinite(temp1)]
        temp1_fit = temp1[np.isfinite(temp1)]
        p0 = np.array([data.max()-data.min(), 50., np.argmax(data), np.nanmin(data)]).astype(float)
        coeff, var_matrix = curve_fit(_gaus, bins_fit, temp1_fit, p0=p0, method='lm', maxfev=50, ftol=1e-5)
        if clipped_wave:
            p0[1] = coeff[1]
            bins_fit = bins[np.isfinite(temp1)]
            temp1_fit = temp1[np.isfinite(temp1)]
            coeff, var_matrix = curve_fit(_gaus, bins_fit, temp1_fit, p0=p0, method='lm', maxfev=50, ftol=1e-5)

        fit_data = _gaus(bins, *coeff)
        chi2 = chisquare(np.array(data, dtype='float64'), f_exp=np.array(fit_data, dtype='float64'))
        fitamplitude = coeff[0]
        fitpos = coeff[1]
        width = coeff[2]*(2.35482)
        base = coeff[3]
        if not (fitamplitude > 1 and fitpos > 0 and width < len(data) and width < fitamplitude and width > 0):
            fitamplitude = np.nan
            fitpos = np.nan
            width = np.nan
            chi2 = np.nan
            base = np.nan
    except RuntimeError:
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
