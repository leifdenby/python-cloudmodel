from pycfd.reference.atmospheric_flow import stratification_profiles
import pysolver.plotting
import common

import numpy as np

profiles = {
    "Isentropic": stratification_profiles.getStandardIsentropicAtmosphere(),
    "Isothermal": stratification_profiles.getStandardIsothermalAtmosphere(),
    "Constant Density": stratification_profiles.getConstantDensityAtmosphere(),
}


plottingHelper = pysolver.plotting.PlottingHelper()
plot = plottingHelper.getPlotter()
plot.suptitle("Hydrostatically balanced atmospheric profiles")

for title, profile in list(profiles.items()):
    z = np.linspace(0.0, 5000.0, 100)

    s = plot.subplot(221)
    plot.xlabel("density/(kg/m3)")
    plot.ylabel("height/m")
    s.plot(profile.rho([z]), z, label=title)
    plot.legend()
    plot.grid(True)
    s = plot.subplot(222)
    plot.ylabel("height/m")
    plot.xlabel("temperature/K")
    s.plot(profile.temp([z]), z, label=title)
    plot.legend()
    plot.grid(True)
    s = plot.subplot(223)
    plot.ylabel("height/m")
    plot.xlabel("potential temperature/K")
    s.plot(profile.pot_temperature([z]), z, label=title)
    plot.legend()
    plot.grid(True)
    s = plot.subplot(224)
    plot.ylabel("height/m")
    plot.xlabel("pressure/(N/m2)")
    s.plot(profile.p([z]), z, label=title)
    plot.legend()
    plot.grid(True)


plot.draw()
input()
