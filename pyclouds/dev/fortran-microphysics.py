# coding: utf-8

import matplotlib.pyplot as plot
import numpy as np
import odespy

from pyclouds import cloud_microphysics, parameterisations
from pyclouds.common import Var, MySolver
from pyclouds.plotting import plot_hydrometeor_evolution

import unified_microphysics.tests.test_common


initial_condition = np.zeros((Var.NUM))
initial_condition[Var.q_v] = 1.35e-2
initial_condition[Var.T] = 297.0
p0 = 99835.0

t_ = np.linspace(0.0, 10.0, 100)

SolverClass = odespy.Fehlberg

constants = unified_microphysics.tests.test_common.um_constants
sol2 = cloud_microphysics.FiniteCondensationTimeMicrophysics(
    r_crit=5e-6, constants=constants
).integrate(initial_condition=initial_condition, t=t_, p0=p0, SolverClass=SolverClass)
sol3 = cloud_microphysics.FortranNoIceMicrophysics(constants=constants).integrate(
    initial_condition=initial_condition, t=t_, p0=p0, SolverClass=SolverClass
)

# sol2.plot()
plot = plot_hydrometeor_evolution(
    [sol2, sol3],
    variables=[
        "q_v",
        "T",
        "q_l",
        "q_r",
        "r_c",
    ],
)
plot.savefig("/scratch/local1/plots/microphysics.png")
