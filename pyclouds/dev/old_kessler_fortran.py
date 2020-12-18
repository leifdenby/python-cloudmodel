# coding: utf-8

import matplotlib.pyplot as plot
import numpy as np
import odespy
import sys

from pyclouds import cloud_microphysics, parameterisations
from pyclouds.common import Var, make_related_constants
from pyclouds.plotting import plot_hydrometeor_evolution

from unified_microphysics.tests.test_common import um_constants


initial_condition = np.zeros((Var.NUM))
q_v = 1.5e-2
q_l = 0.0
T = 288.0

initial_condition[Var.q_v] = q_v
initial_condition[Var.q_l] = q_l

initial_condition[Var.T] = T
p0 = 88676.0  # [Pa]

t_ = [
    0.0,
    0.5,
    1.0,
]  # np.linspace(0., 17., 100)

constants = make_related_constants(um_constants)

q_d = 1.0 - q_v - q_l
q_g = q_d + q_v
Rv, Rd = constants.get("R_v"), constants.get("R_d")
Rg = (Rd * q_d + Rv * q_v) / q_g
print("q_v q_d q_g", q_v, q_d, q_g)
print("pfak = Rg/(Rv * p)*q_g = ", Rg / (Rv * p0) * q_g)
print("qsat_v", parameterisations.pv_sat(T=T) * Rg / (Rv * p0) * q_g)

solutions = []
solutions.append(
    cloud_microphysics.OldATHAMKesslerFortran().integrate(
        initial_condition=initial_condition, t=t_, p0=p0
    )
)


# In[ ]:
