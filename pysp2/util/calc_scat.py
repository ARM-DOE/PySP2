import miepython
import numpy as np
import pandas as pd

from ..data import CalibrationTable

def calc_scattering_psl(size, wavelength=1064.0):
    """
    Calculate the scattering cross section of a PSL sphere of a given a particle size.

    Parameters
    ----------
    size: float
        Diameter of particle in nm.
    wavelength: float
        The wavelength of the laser in nm. Default is 1064 nm (SP2 laser).

    Returns
    -------
    qsca: float
        Scattering efficiency in square meters.
    """

    m = 1.5717
    x = np.pi * size/wavelength
    qext, qsca, qback, g = miepython.mie(m, x)

    return np.pi/4 * (size*1e-9)**2 * qsca

def calc_scattering_psl_table(size):
    """

    Parameters
    ----------
    Calculate the scattering cross section of a PSL sphere of a given a particle size using the lookup Mie table.

    Parameters
    ----------
    size: float
        Diameter of particle in nm.

    Returns
    -------
    qsca: float
        Scattering efficiency in square meters.
    """

    sizes = CalibrationTable["microns"] * 1e3
    scattering = np.zeros_like(size)
    i = 0
    for s in size:
        which_ind = np.argmin(np.abs(sizes - s))
        scattering[i] = CalibrationTable["PSL_1.59-0.0i"][which_ind]
        i+= 1
    return scattering

def calc_scattering_os_table(size):
    """
    Calculate the scattering signature of an organosulfate sphere.

    Parameters
    ----------
    Calculate the scattering cross section of an organosuflate
    sphere of a given a particle size using the lookup Mie table.

    Parameters
    ----------
    size: float
        Diameter of particle in nm.

    Returns
    -------
    qsca: float
        Scattering efficiency in square meters.
    """

    sizes = CalibrationTable["microns"] * 1e3
    which_ind = np.argmin(np.abs(sizes - size))
    return CalibrationTable["Organosulfate_1.48-0.0i"][sizes]

def calc_mass_aquadag(diam):
    """
    Calculate the scattering signature of an organosulfate sphere.

    Parameters
    ----------
    Calculate the scattering cross section of an organosuflate
    sphere of a given a particle size using the lookup Mie table.

    Parameters
    ----------
    size: float
        Diameter of particle in nm.

    Returns
    -------
    mass: float
        Mass in femtograms.
    """
    c = np.array([3.80e-16, -1.99e-13, -5.65e-10, 6.26e-7, -4.82e-5, 0.00422, -0.096])

    return np.polyval(c, diam)