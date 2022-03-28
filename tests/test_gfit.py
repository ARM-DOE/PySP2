import pysp2
import numpy as np


def test_gaussian_fit():
    my_sp2b = pysp2.io.read_sp2(pysp2.testing.EXAMPLE_SP2B)
    my_ini = pysp2.io.read_config(pysp2.testing.EXAMPLE_INI)
    my_binary = pysp2.util.gaussian_fit(my_sp2b, my_ini, parallel=False)
    assert my_binary.PkHt_ch1.max() == 62669.4
    np.testing.assert_almost_equal(
        np.nanmax(my_binary.PkHt_ch0.values), 98708.92915295, decimal=1)
    np.testing.assert_almost_equal(
        np.nanmax(my_binary.PkHt_ch4.values), 54734.05714286, decimal=1)


def test_psds():
    my_sp2b = pysp2.io.read_sp2(pysp2.testing.EXAMPLE_SP2B)
    my_ini = pysp2.io.read_config(pysp2.testing.EXAMPLE_INI)
    my_binary = pysp2.util.gaussian_fit(my_sp2b, my_ini, parallel=False)
    my_hk = pysp2.io.read_hk_file(pysp2.testing.EXAMPLE_HK)
    my_binary = pysp2.util.calc_diams_masses(my_binary)
    ScatRejectKey = my_binary['ScatRejectKey'].values
    assert np.nanmax(
        my_binary['ScatDiaBC50'].values[ScatRejectKey == 0]) < 1000.
    my_psds = pysp2.util.process_psds(my_binary, my_hk, my_ini)
    np.testing.assert_almost_equal(my_psds['NumConcIncan'].max(), 0.95805343)
