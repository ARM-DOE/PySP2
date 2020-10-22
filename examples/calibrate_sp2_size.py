""" This script runs the interactive SP2 calibrator"""
import pysp2
import numpy as np
import pandas as pd

from bokeh.plotting import figure, show, curdoc
from bokeh.layouts import column
from bokeh.models import Button, Range1d
from bokeh.models import Dropdown, ColumnDataSource
from bokeh.server.server import Server

from scipy.optimize import curve_fit

def _gaus(x, a, x0, sigma, base):
    return a * np.exp(-((x - x0)**2/(2 * sigma**2))) + base

if __name__ == "__main__":
    calibration_path = '/nfs/gce/projects/digr/SP2Data4ANL/AquaDag_PSL_cal/Py_generated_numbers/'
    my_cal_data = pysp2.io.read_calibration(calibration_path)
    print(my_cal_data.keys())

    def bkapp(doc):
        def fitButtonCallback(new):
            x_range = [p.x_range.start, p.x_range.end]
            y_range = [p.y_range.start, p.y_range.end]
            x_inds = np.where(np.logical_and(bin_mids >= x_range[0], bin_mids <= x_range[1]), bin_mids, np.nan)
            x = x_inds[np.isfinite(x_inds)]
            y = hist[np.isfinite(x_inds)]
            x_fit = np.linspace(x_range[0], x_range[1], 10000)
            p0 = [y.max(), x[np.argmax(y)], x[np.argmax(y)], 0]
            p_fit = curve_fit(_gaus, x, y, p0, ftol=1e-2)[0]
            y_fit = _gaus(x_fit, *p_fit)
            # p.line(x_fit, y_fit, legend_label='Amp. = %5.2f, Mean = %5.2f, Std. = %5.2f, Base = %5.2f' % p_fit)
            fitDataSource.patch(dict(x=[(slice(0, len(x_fit)), x_fit)], y=[(slice(0, len(y_fit)), y_fit)]))
            fit_str = 'Amp. = %5.2f, Mean = %5.2f, Std. = %5.2f, Base = %5.2f' % (p_fit[0], p_fit[1], p_fit[2], p_fit[3])
            #p.line(x=x_fit, y=y_fit, legend_label=new_data['legend_label'], color='black')
            p.x_range = Range1d(bins[0], bins[-1])
            p.y_range = Range1d(0, hists[0].max())
            p.legend.items = [(settings['sz_key'], [hist_line]), (fit_str, [fit_line])]

        def SzDropDownCallback(new):
            settings['sz_key'] = new.item
            bins = np.linspace(0, 130000, 1000)
            size_datasets = pysp2.util.calc_diams_masses(my_cal_data[settings['sz_key']])
            hist, bins = np.histogram(size_datasets[settings['variable']].where(size_datasets.ScatRejectKey == 0).values,
                                      bins=bins)
            histDataSource.patch(dict(x=[(slice(0, len(bin_mids)), bin_mids)], y=[(slice(0, len(hist)), hist)]))
            p.legend.items = [(settings['sz_key'], [hist_line]), ("", [fit_line])]

        def VariableCallback(new):
            settings['variable'] = new.item
            bins = np.linspace(0, 130000, 1000)
            size_datasets = pysp2.util.calc_diams_masses(my_cal_data[settings['sz_key']])
            hist, bins = np.histogram(size_datasets[settings['variable']].where(size_datasets.ScatRejectKey == 0).values,
                                      bins=bins)
            histDataSource.patch(dict(x=[(slice(0, len(bin_mids)), bin_mids)], y=[(slice(0, len(hist)), hist)]))
            p.legend.items = [(settings['sz_key'], [hist_line]), ("", [fit_line])]

        settings = {'sz_key':'scat_200', 'variable': 'PkHt_ch0'}

        sizes = np.array([350, 300, 269, 220, 200, 170])
        scat_mean = np.ones_like(sizes, dtype=float)

        size_datasets = {}
        i = 0
        # bins = np.logspace(-13, -10, 120)
        bins = np.linspace(0, 130000, 1000)
        hists = np.zeros((len(sizes), len(bins) - 1))
        size_datasets = pysp2.util.calc_diams_masses(my_cal_data[settings['sz_key']])
        hist, bins = np.histogram(size_datasets[settings['variable']].where(size_datasets.ScatRejectKey == 0).values,
                                  bins=bins)
        p = figure(tools="pan,box_zoom,reset,save", x_range=(bins[0], bins[-1]), y_range=(0, hist.max()))
        bin_mids = (bins[:-1] + bins[1:])/2.0
        histDataSource = ColumnDataSource(data=dict(x=bin_mids, y=hist))
        x_fit = np.linspace(0, 130000, 10000)
        fitDataSource = ColumnDataSource(data=dict(x=x_fit, y=np.nan*np.ones_like(x_fit)))
        hist_line = p.line('x', 'y', source=histDataSource, legend_label=settings['sz_key'])
        fit_line = p.line('x', 'y', source=fitDataSource, legend_label="")
        button = Button(label="Curve fit")
        button.on_click(fitButtonCallback)
        size_menu = []
        for keys in my_cal_data.keys():
            size_menu.append((keys, keys))
        SzDropDown = Dropdown(label="Calibration data", menu=size_menu)
        SzDropDown.on_click(SzDropDownCallback)
        vars = ["PkHt_ch0", "PkHt_ch1", "PkHt_ch4", "PkHt_ch5", "FtAmp_ch0", 'FtAmp_ch4']
        vars = [(x, x) for x in vars]
        VarDropDown = Dropdown(label="Variable", menu=vars)
        VarDropDown.on_click(VariableCallback)
        doc.add_root(column(button, SzDropDown, VarDropDown, p))


    server = Server({'/': bkapp}, num_procs=1)
    server.start()
    server.io_loop.add_callback(server.show, "/")
    server.io_loop.start()
