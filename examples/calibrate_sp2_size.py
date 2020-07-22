""" This script runs the interactive SP2 calibrator"""
import pysp2
import numpy as np
from bokeh.plotting import figure, show, curdoc
from bokeh.layouts import column
from bokeh.models import Button

from scipy.optimize import curve_fit

calibration_path = 'C:\\Users\\rjackson\\Documents\\SP2Data4ANL\\AquaDag_PSL_cal\\pyprocessing\\'
my_cal_data = pysp2.io.read_calibration(calibration_path)
calibration_path_igor = 'C:\\Users\\rjackson\\Documents\\SP2Data4ANL\\AquaDag_PSL_cal\\PAPI_Generated_Numbers\\'
my_cal_data_igor = pysp2.io.read_calibration(calibration_path_igor)

size_keys = ['scat_350']
variable = 'PkHt_ch0'
sizes = np.array([350, 300, 269, 220, 200, 170])
scat_mean = np.ones_like(sizes, dtype=float)

size_datasets = {}
i = 0
#bins = np.logspace(-13, -10, 120)
bins = np.linspace(0, 65535, 5000)
hists = np.zeros((len(sizes), len(bins)-1))
for key in size_keys:
    size_datasets[key] = pysp2.util.calc_diams_masses(my_cal_data[key])
    scat_mean[i] = np.nanmedian(size_datasets[key][variable].where(size_datasets[key].ScatRejectKey == 0).values)
    hist, bins = np.histogram(size_datasets[key][variable].where(size_datasets[key].ScatRejectKey == 0).values,
                              bins=bins)
    hists[i] = hist
    i += 1

def _gaus(x, a, x0, sigma, base):
    return a * np.exp(-((x - x0)**2/(2 * sigma**2))) + base

def fitButtonCallback(new):
    x_range = p.x_range
    x_inds = np.where(np.logical_and(bin_mids >= x_range[0], bin_mids <= x_range[1]), bin_mids, np.nan)
    x = x_inds[np.isfinite(x_inds)]
    y = hists[0, np.isfinite(x_inds)]
    x_fit = np.linspace(x_range[0], x_range[1], 10000)
    p0 = [y.max(), x[np.argmax(y)], x[np.argmax(y)], 0]
    p_fit = curve_fit(_gaus, x, y, p0)
    y_fit = _gaus(x_fit, *p_fit)
    #p.line(x_fit, y_fit, legend_label='Amp. = %5.2f, Mean = %5.2f, Std. = %5.2f, Base = %5.2f' % p_fit)
    new_data = dict()
    new_data['x'] = x_fit
    new_data['y'] = y_fit
    new_data['color'] = 'black'
    new_data['legend_label'] = 'Amp. = %5.2f, Mean = %5.2f, Std. = %5.2f, Base = %5.2f' % \
                               (p_fit[0], p_fit[1], p_fit[2], p_fit[3])
    ds.data = new_data
    show(p)


p = figure(tools= "pan,box_zoom,reset,save")
bin_mids = (bins[:-1] + bins[1:])/2.0
p.line(bin_mids, hists[0], legend_label=size_keys[0])
curve = p.line(x=[], y=[])
ds = curve.data_source
button = Button(label="Curve fit")
button.on_click(fitButtonCallback)

show(column(button, p))