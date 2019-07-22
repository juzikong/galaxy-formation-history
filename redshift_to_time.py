def cosmology_calculator(z, H0=69.6, WM=0.286, WV=0.714):
    '''
    http://www.astro.ucla.edu/~wright/CosmoCalc.html
    http://www.astro.ucla.edu/~wright/CC.python
    # symbols   # meaning
    # WR        # Omega(radiation)
    # WK        # Omega curvaturve = 1-Omega(total)
    # age       # age of Universe in units of 1/H0
    # age_Gyr   # value of age in Gyr
    # zage      # age of Universe at redshift z in units of 1/H0
    # zage_Gyr  # value of zage in Gyr
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
