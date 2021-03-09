""" This stores the global variables used by the processing software"""
import numpy as np

class DMTGlobals(object):
    def __init__(self):
        self.ScatMaxPeakHt1 = 60000
        self.ScatMinPeakHt1 = 2500
        self.ScatMaxPeakHt2 = 65535
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
        #mass(fg) = c0Mass2 + c1Mass2 * PkHt_ch5 + c2Mass2 * PkHt_ch5 ^ 2 - Low - gain
        self.c0Mass2 = 0
        self.c1Mass2 = 0.0016815
        self.c2Mass2 = 0
        self.c3Mass2 = 0
        #    // Scattering(1e18
        #    cm2) = c0Scat1 + c1Scat1 * PkHt_ch0 + c2Scat1 * PkHt_ch0 ^ 2 - High - gain
        #    scattering
        #    Variable / g
        self.c0Scat1 = 0
        self.c1Scat1 = 78.141
        self.c2Scat1 = 0
        #    // Scattering(1e18
        #    cm2) = c0Scat2 + c1Scat2 * PkHt_ch4 + c2Scat2 * PkHt_ch4 ^ 2 - Low - gain
        self.c0Scat2 = 0
        self.c1Scat2 = 752.53
        self.c2Scat2 = 0
        self.densitySO4 = 1.8
        self.densityBC = 1.8
        self.TempSTP = 273.15
        self.PressSTP = 1013.25
