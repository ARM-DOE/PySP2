import matplotlib.pyplot as plt
import pysp2
import numpy as np

from scipy.optimize import curve_fit

def linear_fit(x, a):
    return a*x

def pow_law(x, a, b, c):
    return 1000*(a + b*x**c)

Globals = pysp2.util.DMTGlobals()

PkHt_ch4 = np.array([1343.50, 2011.68])
PkHt_ch0 = np.array([20228.93, 13417.36, 3920.25])
Diam0 = np.array([220., 200., 170.])
Diam1 = np.array([200., 220.])

#Scat_not_sat = 1e-18*(Globals.c0Scat1 + Globals.c1Scat1*PkHt_ch0 + Globals.c2Scat1*PkHt_ch0**2)
#Scat_sat = 1e-18*(Globals.c0Scat2 + Globals.c1Scat2*PkHt_ch4 + Globals.c2Scat2*PkHt_ch4**2)

Scat_not_sat = pysp2.util.calc_scattering_psl_table(Diam0)
Scat_sat = pysp2.util.calc_scattering_psl_table(Diam1)
Scat = Scat_sat
Diam = Diam1

# Old fit = 1000*(-0.015256 + 16.835*Scatter**0.15502)

plt.scatter(PkHt_ch4, Scat)
fit = curve_fit(linear_fit, PkHt_ch4, Scat*1e18)
x_range = np.linspace(0, 5000., 100)

plt.plot(x_range, linear_fit(x_range, fit[0])*1e-18, label='y = %5.4ex' % (fit[0]))
plt.plot(x_range, np.polyval([752.53, 0], x_range)*1e-18,
         label='Old calibration')
plt.ylim([0, 5e-12])
plt.legend()
plt.title("Scatter calibration ch4")
plt.xlabel("PkHt_ch4")
plt.ylabel("Scatter [$cm^{-2}$]")
plt.xlim([0, 5000.])
plt.show()
