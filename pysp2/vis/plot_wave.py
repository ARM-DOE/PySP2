import act
import numpy as np
import xarray as xr

from ..util import _gaus


def plot_wave(ds, record_no, chn, plot_fit=True,
              append_to_display=False, init_kwargs=None, **kwargs):
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
    append_to_display : False or ACT HistogramDisplay object
        Default is False. Give an existing ACT HistogramDisplay object to append
        to instead of creating a new one. init_kwargs does not work for an
        existing display object.
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
    time = spectra['time'].values
    inp_data = {}
    inp_data['time'] = xr.DataArray(np.array(time[np.newaxis]),
                                    dims=['time'])
    for ch in range(8):
        inp_data['Data_ch' + str(ch)] = xr.DataArray(
            spectra['Data_ch' + str(ch)].values[np.newaxis, :],
            dims=['time', 'bins'])
    inp_data = xr.Dataset(inp_data)
    bins = np.linspace(0, 100, 100)
    if init_kwargs is None:
        if append_to_display is not False:
            display = append_to_display
        else:
            display = act.plotting.DistributionDisplay(inp_data)
    else:
        display = act.plotting.DistributionDisplay(inp_data, **init_kwargs)
    if 'subplot_index' in kwargs.keys():
        ax = display.plot_size_distribution(
            'Data_ch' + str(chn), bins,
            set_title='Channel %d record %d' % (chn, record_no),
            time=inp_data.time[0],
            subplot_index=kwargs['subplot_index'])
    else:
        ax = display.plot_size_distribution(
            'Data_ch' + str(chn), bins,
            set_title='Channel %d record %d' % (chn, record_no),
            time=inp_data.time[0])

    if plot_fit and chn in [0, 4]:
        xspace = np.linspace(0, 100, 1000)
        amplitude = spectra['FtAmp_ch' + str(chn)].values
        pos = spectra['FtPos_ch' + str(chn)].values
        base = spectra['Base_ch' + str(chn)].values
        width = spectra['PkFWHM_ch' + str(chn)].values
        Y = _gaus(xspace, amplitude, pos, width/2., base)
        ax.plot(xspace, Y)
        ax.text(0.7, 0.5, 'Fit Pos = %3.2f' % pos,
                transform=ax.transAxes)
        ax.text(0.7, 0.6, 'Fit Half Width = %3.2f' % width,
                transform=ax.transAxes)
        ax.text(0.7, 0.7, 'Fit Baseline = %3.2f' % base,
                transform=ax.transAxes)
        ax.text(0.7, 0.8, 'Fit Amplitude = %3.2f' % amplitude,
                transform=ax.transAxes)

    return display


def plot_waves(ds, record_no, plot_fit=True):
    """
    Plots the raw waveforms of all channels for a given record_no.

    Parameters
    ----------
    ds: xarray Dataset
        The dataset to plot the waveform from.
    record_no: int
        The record number to plot.
    plot_fit: bool
        Set to True to plot the gaussian fit. Only used for channels 0 and 4.
    Returns
    -------
    display: ACT HistogramDisplay object
        Returns the ACT
    """

    chns = [i for i in range(8)]
    panel_number = [0, 1, 1, 2, 0, 1, 1, 2]
    legends = [['SCHG ch0', 'ch0_fit', 'SCLG ch4', 'ch4_fit'],
               ['BBHG ch1', 'NBHG ch2', 'BBLG ch5', 'NBLG ch6'],
               ['ch3 HG', 'CH7 LG']]

    if ds['ScatRejectKey'].isel(event_index=record_no) == 3:
        display = plot_wave(ds, record_no, 0, plot_fit=False,
                            init_kwargs={'subplot_shape': (3, )})
        legends[0].remove('ch0_fit')
    else:
        display = plot_wave(ds, record_no, 0, plot_fit=plot_fit,
                            init_kwargs={'subplot_shape': (3, )})
        if plot_fit:
            for t in display.axes[0].texts:
                pos = t.get_position()
                t.set_position((0.1, pos[1]))
    for i, panel in enumerate(panel_number):
        if i != 0:
            display = plot_wave(ds, record_no, chns[i], plot_fit=plot_fit,
                                append_to_display=display,
                                subplot_index=(panel,))
    titles = ['Scattering chanels',
              ' Incandesence channels',
              'Split detector channels']
    xlabels = ['',  '', 'Bin #']
    ylabels = ['Data_ch0, Data_ch4',
               'Data_ch1, Data_ch2\nData_ch5, Data_ch6',
               'Data_ch3, Data_ch7']
    for i, ax in enumerate(display.axes):
        ax.set_title(titles[i])
        ax.set_xlabel(xlabels[i])
        ax.set_ylabel(ylabels[i])
        ax.legend(legends[i])
    return display
