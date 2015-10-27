import numpy as np
import matplotlib.pyplot as plot
from scipy import optimize

from pyclouds import parameterisations
from pyclouds import common

plot.ion()

T_c = np.linspace(-30., 20., 50)
p = 101325.0



pv_sat__atham = parameterisations.ParametersationsWithSpecificConstants(common.ATHAM_constants).pv_sat
pv_sat__ccfm = parameterisations.ParametersationsWithSpecificConstants(common.CCFM_constants).pv_sat
pv_sat__default = parameterisations.ParametersationsWithSpecificConstants(common.default_constants).pv_sat

plot.figure(figsize=(16,30))

plot.subplot(211)
for model in [pv_sat__atham, pv_sat__ccfm, pv_sat__default]:
    plot.plot(T_c, model(T=T_c+273.15), label=str(model), marker='', linestyle=':')
plot.xlabel("Temperature [C]")
plot.ylabel("Saturation vapour pressure [Pa]")
plot.legend(loc='upper left')
plot.grid(True)
plot.title("Saturation vapour pressure with different model constants")

plot.subplot(212)
for model in [pv_sat__atham, pv_sat__ccfm, pv_sat__default]:
    plot.plot(T_c, model.qv_sat(T=T_c+273.15, p=p), label=str(model), marker='', linestyle=':')
plot.xlabel("Temperature [C]")
plot.ylabel("Saturation vapour pressure [Pa]")
plot.legend(loc='upper left')
plot.grid(True)
plot.title("Saturation specific water vapour concentration with different model constants\nAt P=%gPa" % p)

plot.savefig("pv_sat__parameterisation.png")
