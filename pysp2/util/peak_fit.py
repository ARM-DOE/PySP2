import numpy as np
import time
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import chisquare
from .DMTGlobals import DMTGlobals


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

    FtAmp = np.zeros((num_records, 8))
    FtPos = np.zeros((num_records, 8))
    Base = np.zeros((num_records, 8))
    PeakHeight = np.zeros((num_records, 8))
    PeakPos = np.zeros((num_records, 8))
    GaussChiSq = np.zeros((num_records, 8))
    PeakStart = np.zeros((num_records, 8))
    PeakEnd = np.zeros((num_records, 8))
    Width = np.zeros((num_records, 8))
    HalfRise = np.zeros((num_records, 8))
    HalfDecay = np.zeros((num_records, 8))
    Peak2Area = np.zeros((num_records, 8))

    num_trig_pts = int(config['Acquisition']['Pre-Trig Points'])
    start_time = time.time()

    for i in range(num_records):
        if i % 100 == 0 and i != 0:
            print("Records processed: %d / %d" % (i, num_records))
        # Do channel 0 first
        coeffs = _fit_record_gaussian(my_ds, i)
        FtAmp[i, 0] = coeffs['amplitude']
        Base[i, 0] = coeffs['base']
        PeakHeight[i, 0] = coeffs['height']
        PeakPos[i, 0] = coeffs['pos']
        PeakStart[i, 0] = coeffs['start']
        GaussChiSq[i, 0] = coeffs['chi2']
        Width[i, 0] = coeffs['width']

        for chn in [1, 2, 5, 6]:
            coeffs = _fit_record_incan_ave_base(my_ds, i, chn, num_trig_pts)
            Base[i, chn] = coeffs['base']
            PeakHeight[i, chn] = coeffs['height']
            PeakPos[i, chn] = coeffs['pos']
            PeakStart[i, chn] = coeffs['start']
            PeakEnd[i, chn] = coeffs['end']
            HalfRise[i, chn] = coeffs['half_rise']
            HalfDecay[i, chn] = coeffs['half_decay']
            Peak2Area[i, chn] = coeffs['peak2area']

        for chn in [3, 7]:
            coeffs = _split_scatter_fit(my_ds, i, chn)
            Base[i, chn] = coeffs['base']
            PeakHeight[i, chn] = coeffs['height']
            PeakPos[i, chn] = coeffs['pos']
            PeakStart[i, chn] = coeffs['start']

        # Channel 4 is a special snowflake
        coeffs = _gaussian_sat_fit(my_ds, i)
        FtAmp[i, 4] = coeffs['fitamplitude']
        Base[i, 4] = coeffs['base']
        PeakHeight[i, 4] = coeffs['height']
        PeakPos[i, 4] = coeffs['pos']
        FtPos[i, 4] = coeffs['fitpos']
        PeakStart[i, 4] = coeffs['start']
        GaussChiSq[i, 4] = coeffs['chi2']
        Width[i, 4] = coeffs['width']

    # Channel 0
    i = 0
    my_ds['amplitude_ch' + str(i)] = FtAmp[:,i]
    my_ds['amplitude_ch' + str(i)].attrs["long_name"] = "Amplitude for channel %d" % i
    my_ds['base_ch' + str(i)] = Base[:, i]
    my_ds['base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
    my_ds['height_ch' + str(i)] = PeakHeight[:, i]
    my_ds['height_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
    my_ds['width_ch' + str(i)] = Width[:, i]
    my_ds['width_ch' + str(i)].attrs["long_name"] = "Width for channel %d" % i
    my_ds['peak_pos_ch' + str(i)] = PeakPos[:, i]
    my_ds['peak_pos_ch' + str(i)].attrs["long_name"] = "Peak position for channel %d" % i
    my_ds['peak_start_ch' + str(i)] = PeakStart[:, i]
    my_ds['peak_start_ch' + str(i)].attrs["long_name"] = "Peak start for channel %d" % i
    my_ds['chisquare_ch' + str(i)] = GaussChiSq[:, i]
    my_ds['chisquare_ch' + str(i)].attrs["long_name"] = "Chisquare value for channel %d" % i

    # Channels 1, 2, 6, 7
    for i in [1, 2, 5, 6]:
        my_ds['amplitude_ch' + str(i)] = FtAmp[:, i]
        my_ds['amplitude_ch' + str(i)].attrs["long_name"] = "Amplitude for channel %d" % i
        my_ds['base_ch' + str(i)] = Base[:, i]
        my_ds['base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
        my_ds['height_ch' + str(i)] = PeakHeight[:, i]
        my_ds['height_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
        my_ds['half_rise_ch' + str(i)] = HalfRise[:, i]
        my_ds['half_rise_ch' + str(i)].attrs["long_name"] = "Point where rise is at 1/2 height for channel %d" % i
        my_ds['peak2area_ch' + str(i)] = Peak2Area[:, i]
        my_ds['peak2area_ch' + str(i)].attrs["long_name"] = "Peak 2 area for channel %d" % i
        my_ds['half_decay_ch' + str(i)] = HalfDecay[:, i]
        my_ds['half_decay_ch' + str(i)].attrs["long_name"] = "Point where decay is at 1/2 height for channel %d" % i
        my_ds['peak_pos_ch' + str(i)] = PeakPos[:, i]
        my_ds['peak_pos_ch' + str(i)].attrs["long_name"] = "Peak position for channel %d" % i
        my_ds['peak_start_ch' + str(i)] = PeakStart[:, i]
        my_ds['peak_start_ch' + str(i)].attrs["long_name"] = "Peak start for channel %d" % i
        my_ds['peak_end_ch' + str(i)] = PeakEnd[:, i]
        my_ds['peak_end_ch' + str(i)].attrs["long_name"] = "Peak end for channel %d" % i

    for i in [3, 7]:
        my_ds['amplitude_ch' + str(i)] = FtAmp[:, i]
        my_ds['amplitude_ch' + str(i)].attrs["long_name"] = "Amplitude for channel %d" % i
        my_ds['base_ch' + str(i)] = Base[:, i]
        my_ds['base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
        my_ds['height_ch' + str(i)] = PeakHeight[:, i]
        my_ds['height_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
        my_ds['width_ch' + str(i)] = Width[:, i]
        my_ds['width_ch' + str(i)].attrs["long_name"] = "Width for channel %d" % i
        my_ds['peak_pos_ch' + str(i)] = PeakPos[:, i]
        my_ds['peak_pos_ch' + str(i)].attrs["long_name"] = "Peak position for channel %d" % i
        my_ds['peak_split_ch' + str(i)] = PeakStart[:, i]
        my_ds['peak_split_ch' + str(i)].attrs["long_name"] = "Peak split position for channel %d" % i

    # Channel 4
    i = 4
    my_ds['amplitude_ch' + str(i)] = FtAmp[:, i]
    my_ds['amplitude_ch' + str(i)].attrs["long_name"] = "Amplitude for channel %d" % i
    my_ds['base_ch' + str(i)] = Base[:, i]
    my_ds['base_ch' + str(i)].attrs["long_name"] = "Base for channel %d" % i
    my_ds['height_ch' + str(i)] = PeakHeight[:, i]
    my_ds['height_ch' + str(i)].attrs["long_name"] = "Height for channel %d" % i
    my_ds['width_ch' + str(i)] = Width[:, i]
    my_ds['width_ch' + str(i)].attrs["long_name"] = "Width for channel %d" % i
    my_ds['peak_pos_ch' + str(i)] = PeakPos[:, i]
    my_ds['peak_pos_ch' + str(i)].attrs["long_name"] = "Peak position for channel %d" % i
    my_ds['fit_pos_ch' + str(i)] = PeakPos[:, i]
    my_ds['fit_pos_ch' + str(i)].attrs["long_name"] = "Fit position for channel %d" % i
    my_ds['peak_start_ch' + str(i)] = PeakStart[:, i]
    my_ds['peak_start_ch' + str(i)].attrs["long_name"] = "Peak start for channel %d" % i
    my_ds['chisquare_ch' + str(i)] = GaussChiSq[:, i]
    my_ds['chisquare_ch' + str(i)].attrs["long_name"] = "Chisquare value for channel %d" % i
    print(str(num_records) + ' records processed in ' + str(time.time()-start_time) + ' s')
    return my_ds


def _gaus(x, a, x0, sigma, base):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + base


def _fit_record_gaussian(my_ds, record_number):
    """ Only used for channel 0."""
    bins = my_ds.columns.values
    p0 = []
    chn = 0
    data = my_ds['Data_ch' + str(chn)].isel(event_index=record_number).values
    height = np.nan
    pos = np.nan
    start = np.nan
    try:
        coeff, var_matrix = curve_fit(_gaus, bins, data)
        amplitude = coeff[0]
        peakpos = coeff[1]
        width = coeff[2]
        base = coeff[3]
        fit_data = _gaus(bins, *coeff)
        chi2 = chisquare(fit_data, f_exp=data)[1]
    except RuntimeError:
        amplitude = np.nan
        width = np.nan
        peakpos = np.nan
        base = np.nan
        chi2 = np.nan

    if(amplitude < width):
        amplitude = np.nan
        width = np.nan
        peakpos = np.nan
        base = np.nan
        chi2 = np.nan

    if np.isfinite(base):
        peaks, properties = find_peaks(data, width=width)
        try:
            height = properties['prominences'].max()
            pos = peaks[np.argmax(properties['prominences'])]
        except ValueError:
            height = np.nan
            pos = np.nan

        try:
            start = np.where(data - base[chn] < 50)[0][0]
        except IndexError:
            start = np.nan

    # Filter out bad points
    bad = ~np.logical_and.reduce(
        (amplitude > 1, amplitude < 65535, abs(base) < 65535,
         peakpos > 0, peakpos < len(bins), width > 0, width < len(bins)))
    if bad:
        amplitude = np.nan
        base = np.nan
        peakpos = np.nan
        width = np.nan

    fit_coeffs = {'amplitude': amplitude, 'peakpos': peakpos,
                  'width': width, 'base': base, 'chi2': chi2,
                  'height': height, 'pos': pos, 'start': start}
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
    V_maxloc = np.argmax(data)
    peak2area = float(data2.max())/float(data2[20:80].sum())

    if channel in [1, 5]:
        upperbound = 100000.
        # Original code had 0.063, but this seems to be too strict
        lowerbound = 0.063
    else:
        upperbound = 1000000.
        lowerbound = -1000000.


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
        under_half = np.where(diffs < height*0.5)[0]
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
        pos = V_maxloc
    data = data - base
    data = np.abs(data) + data
    if np.isfinite(pos):
        lt5 = np.where(data < 5)[0]
        try:
            start = lt5[lt5 <= pos][-1]
        except IndexError:
            start = np.nan
    fit_coeffs = {'base': base, 'height': height, 'pos': pos, 'start': start,
                  }

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
    error_thrown = False
    start = np.nan
    clipped_wave = False
    global_vars = DMTGlobals()
    data = my_ds['Data_ch' + str(channel)].isel(event_index=record_number).values

    if data.max() - data.min() >= global_vars.ScatMaxPeakHt:
        temp1 = data.astype(float)
        temp1[temp1 == data.max()] = np.nan
        clipped_wave = True

    bins = my_ds.columns.values
    try:
        p0 = [100, 50, 10, 30000]
        coeff, var_matrix = curve_fit(_gaus, bins, data, p0=p0)
        if clipped_wave:
            p0[1] = coeff[1]
            bins_fit = bins[np.isfinite(temp1)]
            temp1_fit = temp1[np.isfinite(temp1)]
            coeff, var_matrix = curve_fit(_gaus, bins_fit, temp1_fit, p0=p0)
        fit_data = _gaus(bins, *coeff)
        chi2 = chisquare(fit_data, f_exp=data)[1]
        fitamplitude = coeff[0]
        fitpos = coeff[1]
        width = coeff[2]
        base = coeff[3]
        if not (fitamplitude > 1 and fitpos > 0 and width < len(data) and width < fitamplitude):
            fitamplitude = np.nan
            fitpos = np.nan
            width = np.nan
            chi2 = np.nan
            base = np.nan
    except RuntimeError as e:
        print(e, record_number)
        error_thrown = True

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
                  'height': height}

    return fit_coeffs





