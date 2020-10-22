import matplotlib.pyplot as plt
import pysp2
import numpy as np

from scipy.optimize import curve_fit

def linear_fit(x, a):
    return a*x


Globals = pysp2.util.DMTGlobals()

PkHt_ch1 = np.array([477.38, 11633.88, 14616.64,])
PkHt_ch5 = np.array([96.41, 1310.64, 2756.61, 5300.])
Diam0 = np.array([69., 179., 201.])
Diam1 = np.array([69., 83., 243., 299.])
Mass0 = pysp2.util.calc_mass_aquadag(Diam0)
Mass1 = pysp2.util.calc_mass_aquadag(Diam1)

#Mass0 = np.pi / 6 * Globals.densityBC * Diam0**3 * 1e15
#Mass1 = np.pi / 6 * Globals.densityBC * Diam1**3 * 1e15


#Scat_not_sat = 1e-18*(Globals.c0Scat1 + Globals.c1Scat1*PkHt_ch0 + Globals.c2Scat1*PkHt_ch0**2)
#Scat_sat = 1e-18*(Globals.c0Scat2 + Globals.c1Scat2*PkHt_ch4 + Globals.c2Scat2*PkHt_ch4**2)

Scat = np.concatenate([PkHt_ch1])
Diam = np.concatenate([Mass0])
# Old fit = 1000*(-0.015256 + 16.835*Scatter**0.15502)
# BC mass (fg) = c0Mass1+c1Mass1*PkHt_ch1+c2Mass1*PkHt_ch1^2 + c3Mass1*PkHt_ch1^3- High-gain incandescent broadband
# BC mass (fg) = c0Mass2+c1Mass2*PkHt_ch5+c2Mass2*PkHt_ch5^2 - Low-gain incandescent broadband

#(sootMass_sat/(0.5236e-9*Globals.densityBC))**(1./3.)
plt.scatter(Scat, Diam)
fit = curve_fit(linear_fit, PkHt_ch1, Diam)

x_range = np.logspace(np.log10(Scat.min()), np.log10(Scat.max()), 100)
#Ch5
#old_coeffs = [0.0016815, 0]
old_coeffs = [0.0001896, 0]
plt.plot(x_range, linear_fit(x_range, fit[0]), label='y = %5.4ex' %
                                                                                      (fit[0]))
plt.plot(x_range, np.polyval(old_coeffs, x_range),
         label='Old calibration')
#plt.plot(x_range*1e18, , label='Old calibration y = 1000*(-0.015256 + 16.835*x**0.15502)')
plt.legend()
plt.title("Black carbon mass ch5 ")
plt.ylabel("Mass [fg]")
plt.xlabel("PkHt_ch5 []")
plt.show()
