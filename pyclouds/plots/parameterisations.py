import numpy as np
import matplotlib.pyplot as plot
from scipy import optimize

from pyclouds import parameterisations

plot.ion()

# from Houghton (1985) via Rogers & Yau (1989)
data_tabulated = np.rec.fromarrays(np.array([
    [-40.0, 1.152e-5,   2.07e-2,    1.62e-5],
    [-30.0, 1.564e-5,   2.16e-2,    1.76e-5],
    [-20.0, 1.616e-5,   2.24e-2,    1.91e-5],
    [-10.0, 1.667e-5,   2.32e-2,    2.06e-5],
    [0.0,   1.717e-5,   2.40e-2,    2.21e-5],
    [10.,   1.766e-5,   2.48e-2,    2.36e-5],
    [20.,   1.815e-5,   2.55e-2,    2.52e-5],
    [30.,   1.862e-5,   2.63e-2,    2.69e-5],
]).T,
dtype=[
    ('temp', 'f8'), # temperature [K]
    ('mu', 'f8'), # dynamic viscosity [kg m-1 s-1]
    ('K', 'f8'), # thermal conductivity of air [J m-1 s-1 K-1]
    ('D', 'f8')] # coefficient of diffusion of water vapour in air [m2 s-1]
)


dT = np.max(data_tabulated.temp) - np.min(data_tabulated.temp)
T_c = np.linspace(np.min(data_tabulated.temp) - 0.2*dT, np.max(data_tabulated.temp) + 0.2*dT)

p = 101325.0


plot.figure(figsize=(10,20))

Ka = parameterisations.Ka

plot.subplot(211)
plot.plot(T_c, Ka(T_c+273.15), label=str(Ka), marker='', linestyle=':')
plot.plot(data_tabulated.temp, data_tabulated.K, label="Rogers & Yau", marker='s', linestyle='')
plot.xlabel("Temperature [C]")
plot.ylabel("Thermal conductivity of air [J m-1 s-1 K-1]")
plot.legend(loc='lower right')
plot.grid(True)
plot.title("Thermal conductivity of air")

plot.subplot(212)

Dv_ATHAM_constants = parameterisations.WaterVapourDiffusionCoefficient.ATHAM_constants
Dv_ATHAM = parameterisations.WaterVapourDiffusionCoefficient(constants=Dv_ATHAM_constants)

# assumed powerlaw, tag log so that we can use polynomial fitting from numpy
b_, a_ = np.polyfit(np.log((data_tabulated.temp+273.15)/273.15), np.log(data_tabulated.D), 1)
b = b_
a = np.exp(a_)
Dv__fitted = parameterisations.WaterVapourDiffusionCoefficient(constants={'a': a, 'b': b})
Dv = parameterisations.WaterVapourDiffusionCoefficient()

plot.grid(True)

for model in [Dv_ATHAM, Dv__fitted, Dv]:
    plot.plot(T_c, model(T_c+273.15, p=p), label=str(model), marker='', linestyle=':')
plot.plot(data_tabulated.temp, data_tabulated.D, label="Rogers & Yau", marker='s', linestyle='')
plot.xlabel("Temperature [C]")
plot.ylabel("Coefficient of diffusion of water vapour in air [m2 s-1]")
plot.legend(loc='lower right')
plot.grid(True)
plot.title("Diffusion of water vapour in air (at %g Pa)" % p)

plot.savefig("parameterisations.png")