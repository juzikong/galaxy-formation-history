import matplotlib.pyplot as plt
import numpy as np
import math


def cosmology_calculator(z, H0=69.6, WM=0.286, WV=0.714):
    '''
    http://www.astro.ucla.edu/~wright/CosmoCalc.html
    http://www.astro.ucla.edu/~wright/CC.python
    # symbols   # meaning
    # H0        # Hubble constant
    # WM        # Omega_M
    # WV        # Omega_vac
    # z         # redshift
    # WR        # Omega(radiation)
    # WK        # Omega curvaturve = 1-Omega(total)
    # age       # age of Universe in units of 1/H0
    # age_Gyr   # value of age in Gyr
    # zage      # age of Universe at redshift z in units of 1/H0
    # zage_Gyr  # value of zage in Gyrs
    # a         # 1/(1+z), the scale factor of the Universe
    # az        # 1/(1+z(object))
    # Tyr       # coefficent for converting 1/H into Gyr
    '''
    Tyr = 977.8  # coefficent for converting 1/H into Gyr
    h = H0 / 100.
    WR = 4.165E-5 / (h * h)  # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1 - WM - WR - WV
    az = 1.0 / (1 + 1.0 * z)
    age = 0.
    n = 1000  # number of points in integrals
    for i in range(n):
        a = az * (i + 0.5) / n
        adot = (WK + (WM / a) + (WR / (a * a)) + (WV * a * a)) ** 0.5
        age = age + 1. / adot
    zage = az * age / n
    zage_Gyr = (Tyr / H0) * zage
    return zage_Gyr


H0 = 69.6
WM = 0.286
WV = 0.714
print("Standard age of the universe:", cosmology_calculator(0, H0, WM, WV))

# ### plot redshift--time relation ###
# redshift_list = np.arange(0.003, 0.1, 0.001)
# redshift_list = np.append(redshift_list, np.arange(0.1, 10, 0.01))
# redshift_list = np.append(redshift_list, np.arange(10, 100, 0.1))
# time_list = cosmology_calculator(redshift_list, H0, WM, WV)
# plt.yscale("log")
# plt.xlim(0, math.ceil(cosmology_calculator(0, H0, WM, WV)))
# plt.ylim(0.001, 100)
# plt.xlabel("time [Gyr]")
# plt.ylabel("redshift")
# plt.title("age of the universe: {} Gyr".format(round(cosmology_calculator(0), 2)))
# plt.plot(time_list, redshift_list)
# plt.show()
