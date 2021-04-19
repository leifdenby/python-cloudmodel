# coding: utf-8

import matplotlib.pyplot as plot
import numpy as np
import odespy

from pyclouds import cloud_microphysics, parameterisations
from pyclouds.common import Var, MySolver
from pyclouds.plotting import plot_hydrometeor_evolution


initial_condition = np.zeros((Var.NUM))
initial_condition[Var.q_v] = 0.017
initial_condition[Var.T] = 285.0
p0 = 101325.1  # [Pa]

t_ = np.linspace(0.0, 10.0, 100)

SolverClass = odespy.Fehlberg

microphysics_model = cloud_microphysics.MoistAdjustmentMicrophysics()
sol = microphysics_model.integrate(initial_condition=initial_condition, t=t_, p0=p0)
# sol1 = microphysics_model.integrate(initial_condition=initial_condition, t=t_, p0=p0, iterations=3)
# sol2 = cloud_microphysics.FiniteCondensationTimeMicrophysics().integrate(initial_condition=initial_condition, t=t_, p0=p0, SolverClass=SolverClass)
sol2 = cloud_microphysics.FiniteCondensationTimeMicrophysics(r_crit=6e-6).integrate(
    initial_condition=initial_condition, t=t_, p0=p0, SolverClass=SolverClass
)
sol3 = cloud_microphysics.FiniteCondensationTimeMicrophysics(r_crit=np.inf).integrate(
    initial_condition=initial_condition, t=t_, p0=p0, SolverClass=SolverClass
)

# sol2.plot()
plot = plot_hydrometeor_evolution([sol2, sol3], variables=["q_v", "r_c", "q_l", "q_r"])
plot.savefig("/scratch/local1/plots/microphysics.png")
