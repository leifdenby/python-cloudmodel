"""
Studying the importance of the temperature dependence of the moist lapse rate.

Equations and constants from http://en.wikipedia.org/wiki/Lapse_rate

NB: The expressions below are not giving sensible answers - 24/2/2014
"""

import numpy as np
from matplotlib import pyplot as plot


g = 9.8076  # m/s^2
H_v = 2260000.0  # J/g
r = 0.6219897  # g/g
R = 8.314  # J/mol/K

R_sd = 287.0  # J/g/K
R_sw = 462.0  # J/g/K
eps = R_sd / R_sw
c_pd = 1003.5  # J/g/K


def lapse_wet(T):
    return (
        g
        * (1.0 + (H_v * r) / (R_sd * T))
        / (c_pd + H_v ** 2.0 * r / (R_sw * T ** 2.0))
        * 1000.0
    )


print(lapse_wet(273.0 + 40.0), lapse_wet(273.0), lapse_wet(273.0 - 40.0))
