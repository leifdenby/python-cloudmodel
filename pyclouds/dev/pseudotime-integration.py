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

output_filename = "/tmp/cloud-profile-time-integration.png"

ambient_stratification = stratification_profiles.Soong1973Dry()
T_e = lambda z: ambient_stratification.temp(z)
p_e = ambient_stratification.p


environment = ambient_stratification

t = numpy.linspace(0, 1200.0, 800)

w0 = 0.1
T0 = T_e(0.0) + 0.2
z0 = 50.0
r0 = 1000.0
initial_condition = Var.make_state(r=r0, w=w0, T=T0, q_l=0.001, q_v=0.015, z=z0)


stopping_criterion = None

beta = 0.0
D = 0.0
gamma = 1.0


profiles = []

microphysics = cloud_microphysics.FiniteCondensationTimeMicrophysics(disable_rain=False)
cloud_model = cloud_equations.FullThermodynamicsCloudEquations(
    gamma=gamma, D=D, beta=beta, environment=environment, microphysics=microphysics
)
profiles.append(cloud_model.integrate_in_time(initial_condition, t))

microphysics = cloud_microphysics.FiniteCondensationTimeMicrophysics(disable_rain=True)
cloud_model = cloud_equations.FullThermodynamicsCloudEquations(
    gamma=gamma, D=D, beta=beta, environment=environment, microphysics=microphysics
)
profiles.append(cloud_model.integrate_in_time(initial_condition, t))

fig = plotting.plot_profiles(
    profiles,
    variables=["r", "w", "T", "q_v", "q_l", "q_r"],
    initial_condition=initial_condition,
    labels_ncol=1,
)
fig.savefig(output_filename)
print("output written to %s" % output_filename)
