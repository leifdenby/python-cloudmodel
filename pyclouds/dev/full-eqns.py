# coding: utf-8
import matplotlib

matplotlib.use("Agg")

import numpy as np
import odespy, numpy

from pycfd.reference.atmospheric_flow import (
    stratification_profiles,
    saturation_calculation,
)
from pyclouds import cloud_equations, plotting, cloud_microphysics
from pyclouds.common import Var

output_filename = "/tmp/cloud-profile-comparison.png"

ambient_stratification = stratification_profiles.Soong1973Dry()
T_e = lambda z: ambient_stratification.temp(z)
p_e = ambient_stratification.p


environment = ambient_stratification
z_max = 5e3
z_points = numpy.linspace(100.0, z_max, 100)

w0 = 0.1
T0 = T_e(0.0) + 0.2
# p, r, w, T, q_v, q_r, q_l, q_i
initial_condition = Var.make_state(r=500, w=w0, T=T0, q_l=0.001, q_v=0.015)


stopping_criterion = None

beta = 0.0
d_beta = 0.02
D = 0.0
dD = 0.0
gamma = 1.0

microphysics = cloud_microphysics.FiniteCondensationTimeMicrophysics(disable_rain=True)

default_kwargs = {
    "environment": environment,
    "microphysics": microphysics,
}

integration_kwargs = {
    "initial_condition": initial_condition,
    "z": z_points,
}

profiles = []

cloud_model = cloud_equations.FullThermodynamicsCloudEquations(
    gamma=gamma, D=D + dD, beta=beta + d_beta, **default_kwargs
)
profiles.append(cloud_model.integrate(initial_condition, z_points))

cloud_model = cloud_equations.FullEquationsSatMicrophysics(
    gamma=gamma, D=D + dD, beta=beta + d_beta, **default_kwargs
)
profiles.append(cloud_model.integrate(initial_condition, z_points))

cloud_model = cloud_equations.DryAirOnly(
    gamma=gamma, D=D + dD, beta=0.0, **default_kwargs
)
profiles.append(cloud_model.integrate(initial_condition, z_points))

fig = plotting.plot_profiles(
    profiles,
    variables=["r", "w", "T", "q_v", "T__tephigram", "q_l"],
    initial_condition=initial_condition,
)
fig.savefig(output_filename)
print("output written to %s" % output_filename)
