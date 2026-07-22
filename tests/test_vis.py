import matplotlib
import pytest
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

import pysp2
from pysp2.util.normalized_derivative_method import plot_normalized_derivative
from pysp2.util.normalized_derivative_method import MLEConfig, mle_tau_moteki_kondo, compute_d2_moteki_kondo
from pysp2.util.normalized_derivative_method import compute_sigma_moteki_kondo
from pysp2.util.normalized_derivative_method import plot_incident_irradiance
from pysp2.util.normalized_derivative_method import plot_scattering_cross_section
from pysp2.vis.plot_wave import plot_wave

matplotlib.use("Agg")

my_sp2b = pysp2.io.read_sp2(pysp2.testing.EXAMPLE_SP2B_PSL, arm_convention=False)
my_ini = pysp2.io.read_config(pysp2.testing.EXAMPLE_INI_PSL)
event = 152

@pytest.mark.mpl_image_compare(tolerance=10)
def test_plot_normalized_derivative():
    my_binary = pysp2.util.gaussian_fit(my_sp2b, my_ini, parallel=False, baseline_to_zero=True)
    dSdt_norm = pysp2.util.central_difference(my_binary, normalize=True, baseline_to_zero=True)

    # Test the plotting function for channel 0 and record number 720
    ax = plot_normalized_derivative(my_binary, dSdt_norm, record_no=event, chn=0, plot_scattering_signal=True)
    fig = ax.figure
    
    return fig

@pytest.mark.mpl_image_compare(tolerance=10)
def test_plot_wave():

    my_binary = pysp2.util.gaussian_fit(my_sp2b, my_ini, parallel=False, baseline_to_zero=False)

    # Test the plotting function for channel 0 and record number 720
    display = plot_wave(my_binary, record_no=event, chn=0, plot_fit=True)
    fig = display.axes[0].figure
    return fig

@pytest.mark.mpl_image_compare(tolerance=10)
def test_plot_incident_irradiance():

    my_binary = pysp2.util.gaussian_fit(my_sp2b, my_ini, parallel=False, baseline_to_zero=True)
    dSdt = pysp2.util.central_difference(my_binary, normalize=True, baseline_to_zero=True)
    print("dSdt dimensions:", dSdt.dims)

    cfg = MLEConfig(
    h=0.4,           # example: 0.4 microseconds
    sigma_bar= 18.5*0.4,  # example; use your measured average width
    delta_sigma=1.2*0.4, # example; use your measured width std dev
    A1=0.37*2.44,
    A2=(1.6e-2)*2.44**(1/2),
    A3=6.2e-4,
)
    
    tau = mle_tau_moteki_kondo(
        S=my_binary,
        norm_deriv=dSdt,
        p=13,
        ch="Data_ch0",
        event_index=event,
        min_start=15,
        width_metric="fwtm",
        config=cfg,
    )

    d2 = compute_d2_moteki_kondo(
        S=my_binary,
        norm_deriv=dSdt,
        tau_hat=tau,
        p=13,
        ch="Data_ch0",
        event_index=event,
        min_start=15,
        width_metric="fwtm",
        config=cfg,
    )

    sigma_ds = compute_sigma_moteki_kondo(
        S=my_binary,
        norm_deriv=dSdt,
        tau_hat=tau,
        d2=d2,
        p=13,
        ch="Data_ch0",
        event_index=event,
        min_start=15,
        width_metric="fwtm",
        config=cfg,
    )

    # Test the plotting function for channel 0 and record number 720
    ax = plot_incident_irradiance(
        S=my_binary,
        ds=dSdt,
        record_no=event,
        chn=0,
        plot_scattering_signal=True,
        sigma_ds=sigma_ds,
        time_units="us",
    )

    fig = ax.figure
    
    return fig

@pytest.mark.mpl_image_compare(tolerance=10)
def test_plot_scattering_cross_section():

    my_binary = pysp2.util.gaussian_fit(my_sp2b, my_ini, parallel=False, baseline_to_zero=True)
    dSdt = pysp2.util.central_difference(my_binary, normalize=True, baseline_to_zero=True)
    print("dSdt dimensions:", dSdt.dims)

    cfg = MLEConfig(
    h=0.4,           # example: 0.4 microseconds
    sigma_bar= 18.5*0.4,  # example; use your measured average width
    delta_sigma=1.2*0.4, # example; use your measured width std dev
    A1=0.37*2.44,
    A2=(1.6e-2)*2.44**(1/2),
    A3=6.2e-4,
)
    
    tau = mle_tau_moteki_kondo(
        S=my_binary,
        norm_deriv=dSdt,
        p=13,
        ch="Data_ch0",
        event_index=event,
        min_start=15,
        width_metric="fwtm",
        config=cfg,
    )

    d2 = compute_d2_moteki_kondo(
        S=my_binary,
        norm_deriv=dSdt,
        tau_hat=tau,
        p=13,
        ch="Data_ch0",
        event_index=event,
        min_start=15,
        width_metric="fwtm",
        config=cfg,
    )

    sigma_ds = compute_sigma_moteki_kondo(
        S=my_binary,
        norm_deriv=dSdt,
        tau_hat=tau,
        d2=d2,
        p=13,
        ch="Data_ch0",
        event_index=event,
        min_start=15,
        width_metric="fwtm",
        config=cfg,
    )

    # Test the plotting function for channel 0 
    fig = plot_scattering_cross_section(
        S=my_binary,
        sigma_ds=sigma_ds,
        record_no=event,
        chn=0,
        ch="Data_ch0",
        h=0.4,
        time_units="us"
        )

    return fig