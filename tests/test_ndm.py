import pysp2
import numpy as np
np.set_printoptions(threshold=np.inf)

from pysp2.util.normalized_derivative_method import MLEConfig, mle_tau_moteki_kondo, compute_d2_moteki_kondo
from pysp2.util.normalized_derivative_method import compute_sigma_moteki_kondo, compute_normalized_incident_irradiance_moteki_kondo
event=152
my_sp2b = pysp2.io.read_sp2(pysp2.testing.EXAMPLE_SP2B_PSL, arm_convention=False)
my_ini = pysp2.io.read_config(pysp2.testing.EXAMPLE_INI_PSL)
my_binary = pysp2.util.gaussian_fit(my_sp2b, my_ini, parallel=False)

def test_central_difference():
    dSdt = pysp2.util.central_difference(my_binary, normalize=False, baseline_to_zero=False)
    
    np.testing.assert_almost_equal(dSdt['Data_ch4'].isel(event_index=event, time=0).item(), 
                                   2.50, decimal=2)
    np.testing.assert_almost_equal(dSdt['Data_ch4'].isel(event_index=event, time=99).item(), 
                                   13.33, decimal=2)
    np.testing.assert_almost_equal(dSdt['Data_ch4'].isel(event_index=event, time=2).item(), 
                                    0.83, decimal=2)
    assert np.isfinite(dSdt).all()

    dSdt_norm = pysp2.util.central_difference(my_binary, normalize=True, baseline_to_zero=True)
    np.testing.assert_almost_equal(dSdt_norm['Data_ch0'].isel(event_index=event, time=0).item(), 
                                   -0.70, decimal=2)
    np.testing.assert_almost_equal(dSdt_norm['Data_ch0'].isel(event_index=event, time=99).item(), 
                                   0.799, decimal=2)
    np.testing.assert_almost_equal(dSdt_norm['Data_ch0'].isel(event_index=event, time=2).item(), 
                                    -0.078, decimal=2)
    
def test_ndm_moteki_kondo():
    dSdt = pysp2.util.central_difference(my_binary, normalize=True, baseline_to_zero=True)  

    cfg = MLEConfig(
    h=0.4,           # example: 0.4 microseconds
    sigma_bar= 18.5*0.4,  # example; use your measured average width
    delta_sigma=1.2*0.4, # example; use your measured width std dev
    A1=0.37*2.44,
    A2=(1.6e-2)*2.44**(1/2),
    A3=6.2e-4,
)
    
    ## Test one event ##################################################
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

    tau_val_true = my_binary['Data_ch0'].isel(event_index=event).argmax().item()*0.4
    # Test that the estimated tau for a subset of results is close to the true value for the event
    for i in range(10, 15):
        np.testing.assert_allclose(tau[i], tau_val_true, atol=0.3)
    
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

    k_min_local = int(d2.argmin(dim="k").item())
    # Get corresponding tau value
    tau_best = tau.isel(k=k_min_local).item()

    # Assert closeness
    np.testing.assert_allclose(
        tau_best,
        tau_val_true,
        atol=0.3,  # absolute tolerance = 0.3 microseconds
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

    # example: use the best sigma value from your analysis, divided by 4.29193 to convert FWTM value of
    # 18.51*np.sqrt(np.log(10)/np.log(2)) to sigma where 33.7366 is the average FWHM in us
    sigma_best = (33.7366*0.4)/4.29193

    np.testing.assert_allclose(
         sigma_ds['sigma_hat'].values,
         sigma_best,
         atol=0.12,  # absolute tolerance = 1.5 microseconds
    )
    
    # Test the normalized irradiance function 
    I_norm = compute_normalized_incident_irradiance_moteki_kondo(
        sigma_out=sigma_ds,
    )

    y_scatter_background_shifted = (
        my_binary['Data_ch0'].isel(event_index=event).values - 
        np.nanmin(my_binary['Data_ch0'].isel(event_index=event).values)
    )

    # test for peak area only
    for i in range(15,75):
        np.testing.assert_allclose(
            (I_norm * np.nanmax(y_scatter_background_shifted))[i],
            y_scatter_background_shifted[i],
            atol=4500,  # absolute tolerance ~ 10% of the max scattering signal value
        )