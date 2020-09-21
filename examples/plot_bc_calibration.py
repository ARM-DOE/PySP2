import matplotlib.pyplot as plt
import pysp2
import numpy as np

from scipy.optimize import curve_fit


Globals = pysp2.util.DMTGlobals()

PkHt_ch1 = np.array([12786.38, 13262.56, 13511.74, 11633.88, 14616.64, ])
PkHt_ch5 = np.array([1188., 1490.35, 2756.61, 5430.])
Diam0 = np.array([83., 102., 130,  179., 201.])*1e-7
Diam1 = np.array([179., 201., 243., 299])*1e-7
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
fit = np.polyfit(Scat, Diam, 2)
fit = np.append([0], fit)
x_range = np.logspace(np.log10(Scat.min()), np.log10(Scat.max()), 100)

plt.plot(x_range, np.polyval(fit, x_range), label='y = + %5.4ex**3 +  %5.4ex**2 + %5.4ex + %5.4e ' %
                                                                                      (fit[0], fit[1], fit[2], fit[3]))
plt.plot(x_range, np.polyval([Globals.c3Mass1, Globals.c2Mass1, Globals.c1Mass1, Globals.c0Mass1], x_range),
         label='Old calibration')
#plt.plot(x_range*1e18, , label='Old calibration y = 1000*(-0.015256 + 16.835*x**0.15502)')
plt.legend()
plt.title("Black carbon mass hi-gain ")
plt.ylabel("Mass [fg]")
plt.xlabel("PkHt_ch1 []")
plt.show()
