""" This stores the global variables used by the processing software"""
import numpy as np
import pandas as pd


class DMTGlobals(object):
    """
    This class stores the default parameters used for particle processing.

    Attributes
    ----------
    ScatMaxPeakHt1: int
        The maximum peak value for scattering channel 0.
    ScatMinPeakHt1: int
        The minimum peak value for scattering channel 0
    ScatMaxPeakHt2: int
        The maximum peak value for scattering channel 4.
    ScatMinPeakHt2: int
        The minimum peak value for scattering channel 4
    ScatMinWidth: int
        The minimum peak width for the scattering channels
    ScatMaxWidth: int
        The maximum peak width for the scattering channels
    IncanMinPeakHt1: int
        The minimum peak height for channel 1
    IncanMinPeakHt2: int
        The minimum peak height for channel 5
    IncanMidWidth: int
        The minimum width for the incadescence channels
    IncanMaxWidth: int
        The maximum width for the incadescence channels
    IncanMinPeakPos: int
        The minimum peak position for the incadescence channels.
    IncanMaxPeakPos: int
        The maximum peak position for the incadescence channels.
    IncanMinPeakRatio: float
        The minimum peak ch5/ch6 peak ratio.
    IncanMaxPeakRatio: float
        The maximum peak ch5/ch6 peak ratio.
    c0Mass1, c1Mass1, c2Mass1: floats
        Calibration mass coefficients for ch1
    c0Mass2, c1Mass2, c2Mass2: floats
        Calibration mass coefficients for ch2
    c0Scat1, c1Scat1, c2Scat1: floats
        Calibration scattering coefficients for ch0
    c0Scat2, c1Scat2, c2Scat2: floats
        Calibration scattering coefficients for ch4
    densitySO4, BC: float
        Density of SO4, BC
    tempSTP, presSTP: float
        Temperature [Kelvin] and pressure [hPa] at STP.
    """

    def __init__(self, cal_file_name=None):
        self.ScatMaxPeakHt1 = 60000
        self.ScatMinPeakHt1 = 250
        self.ScatMaxPeakHt2 = 60000
        self.ScatMinPeakHt2 = 250
        self.ScatMinWidth = 10
        self.ScatMaxWidth = 90
        self.ScatMinPeakPos = 20
        self.ScatMaxPeakPos = 90
        self.IncanMinPeakHt1 = 200
        self.IncanMinPeakHt2 = 200
        self.IncanMaxPeakHt1 = 60000
        self.IncanMaxPeakHt2 = 60000
        self.IncanMinWidth = 5
        self.IncanMaxWidth = np.inf
        self.IncanMinPeakPos = 20
        self.IncanMaxPeakPos = 90
        self.IncanMinPeakRatio = 0.1
        self.IncanMaxPeakRatio = 25
        self.IncanMaxPeakOffset = 11
        # Default values are for Brookhaven SP2, unit  # 24 - 2010-11-17 - RS
        #  mass(fg) = c0Mass1 + c1Mass1 * PkHt_ch1 + c2Mass1 * PkHt_ch1 ^ 2 + c3Mass1 * PkHt_ch1 ^ 3 - High - gain
        self.c0Mass1 = 0
        self.c1Mass1 = 0.0001896
        self.c2Mass1 = 0
        self.c3Mass1 = 0
        self.c0Mass2 = 0
        self.c1Mass2 = 0.0016815
        self.c2Mass2 = 0
        self.c3Mass2 = 0
        self.c0Scat1 = 0
        self.c1Scat1 = 78.141
        self.c2Scat1 = 0
        self.c0Scat2 = 0
        self.c1Scat2 = 752.53
        self.c2Scat2 = 0
        self.densitySO4 = 1.8
        self.densityBC = 1.8
        self.TempSTP = 273.15
        self.PressSTP = 1013.25
        if cal_file_name is not None:
            df = pd.read_csv(cal_file_name, sep='\t')
            for i in range(len(df["CalName"].values)):
                setattr(self, df["CalName"].values[i], df["CalValue"].values[i])
            del df
