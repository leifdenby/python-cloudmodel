import numpy as np

from pyclouds import cloud_microphysics
from pyclouds.common import Var

initial_condition = np.zeros((Var.NUM))
initial_condition[Var.q_v] = 0.03
initial_condition[Var.T] = 300.
p0 = 101325.1  # [Pa]

t_ = np.linspace(0., 10., 100)

microphysics_model = cloud_microphysics.MoistAdjustmentMicrophysics()
sol = microphysics_model.integrate(initial_condition=initial_condition, t=t_, p0=p0)
sol2 = cloud_microphysics.FiniteCondensationTimeMicrophysics().integrate(initial_condition=initial_condition, t=t_, p0=p0)

sol2.plot()
