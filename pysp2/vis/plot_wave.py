import act
import numpy as np

from ..util import _gaus


def plot_wave(ds, record_no, chn, plot_fit=True, init_kwargs=None, **kwargs):
    """
    Plots the raw waveform for a given record_no and channel.

    Parameters
    ----------
    ds: xarray Dataset
        The dataset to plot the waveform from.
    record_no: int
        The record number to plot.
    chn: int
        The channel to plot.
    plot_fit: bool
        Set to True to plot the gaussian fit. Only used for channels 0 and 4.
    init_kwargs: dict
        Dictionary of keyword arguments to pass into initialization of
        :func:`act.plotting.display.HistogramDisplay` object
    kwargs:
        Additional keyword arguments are passed into :func:`matplotlib.pyplot.step`

    Returns
    -------
    display: ACT HistogramDisplay object
        Returns the ACT
    """
    spectra = ds.isel(event_index=record_no)
    bins = np.linspace(0, 100, 100)
    if init_kwargs is None:
        display = act.plotting.HistogramDisplay(spectra)
    else:
        display = act.plotting.HistogramDisplay(spectra, **init_kwargs)
    ax = display.plot_size_distribution('Data_ch' + str(chn), bins, set_title='Channel %d record %d' % (chn, record_no))

    if plot_fit and chn in [0, 4]:
        xspace = np.linspace(0, 100, 1000)
        amplitude = spectra['FtAmp_ch' + str(chn)].values
        pos = spectra['FtPos_ch' + str(chn)].values
        base = spectra['Base_ch' + str(chn)].values
        width = spectra['PkFWHM_ch' + str(chn)].values
        Y = _gaus(xspace, amplitude, pos, width/2., base)
        ax.plot(xspace, Y)
        ax.text(0.7, 0.5, 'Fit Pos = %3.2f' % pos, transform=ax.transAxes)
        ax.text(0.7, 0.6, 'Fit Half Width = %3.2f' % width, transform=ax.transAxes)
        ax.text(0.7, 0.7, 'Fit Baseline = %3.2f' % base, transform=ax.transAxes)
        ax.text(0.7, 0.8, 'Fit Amplitude = %3.2f' % amplitude, transform=ax.transAxes)

    return display
