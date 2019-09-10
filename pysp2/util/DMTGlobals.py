""" This stores the global variables used by the processing software"""
import numpy as np

class DMTGlobals(object):
    def __init__(self):
        self.ScatMaxPeakHt = 60000
        self.ScatMinPeakHt = 2500
        self.ScatMinWidth = 10
        self.ScatMaxWidth = 90
        self.ScatMinPeakPos = 20
        self.ScatMaxPeakPos = 90
        self.IncanMinPeakHt = 20
        self.IncanMaxPeakHt = 60000
        self.IncanMinWidth = 5
        self.IncanMaxWidth = np.inf
        self.IncanMinPeakPos = 20
        self.IncanMaxPeakPos = 90
        self.IncanMinPeakRatio = 0.1
        self.IncanMaxPeakRatio = 25
        self.IncanMaxPeakOffset = 11

