import numpy as np
import matplotlib.pyplot as plot
from scipy import optimize

from ..reference import parameterisations
from .. import Var
from ..reference import constants

savefig = False

plot_kwargs = {}
if savefig:
    plot.ion()
    plot_kwargs["figsize"] = (16, 30)

T_c = np.linspace(-30.0, 20.0, 50)

pv_sat__atham = parameterisations.ParametersationsWithSpecificConstants(
    constants.ATHAM_constants
).pv_sat
pv_sat__ccfm = parameterisations.ParametersationsWithSpecificConstants(
    constants.CCFM_constants
).pv_sat
pv_sat__default = parameterisations.ParametersationsWithSpecificConstants(
    constants.default_constants
).pv_sat

plot.figure(**plot_kwargs)

plot.subplot(211)
for model in [pv_sat__atham, pv_sat__ccfm, pv_sat__default]:
    plot.plot(T_c, model(T=T_c + 273.15), label=str(model), marker="", linestyle=":")
plot.xlabel("Temperature [C]")
plot.ylabel("Saturation vapour pressure [Pa]")
plot.legend(loc="upper left")
plot.grid(True)
plot.title("Saturation vapour pressure with different model constants")

plot.subplot(212)
p = 88676.0  # 101325.0
for model in [pv_sat__atham, pv_sat__ccfm, pv_sat__default]:
    plot.plot(
        T_c,
        model.qv_sat(T=T_c + 273.15, p=p),
        label=str(model),
        marker="",
        linestyle=":",
    )
plot.xlabel("Temperature [C]")
plot.ylabel("Saturation specific water concentration [kg/kg]")
plot.legend(loc="upper left")
plot.title("Ambient pressure p=%gPa" % p)
plot.grid(True)

plot.suptitle(
    "Saturation specific water vapour concentration with different model constants"
)

if savefig:
    plot.savefig("pv_sat__parameterisation.png")
else:
    plot.draw()
    plot.show()
