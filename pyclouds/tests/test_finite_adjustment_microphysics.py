import numpy as np

from pyclouds import cloud_microphysics
from pyclouds.common import Var


def test_integration():
    initial_condition = Var.make_state(T=300., p=101325.0, q_v=0.0)

    t_ = np.linspace(0., 10., 100)

    microphysics_model = cloud_microphysics.MoistAdjustmentMicrophysics()
    sol = microphysics_model.integrate(initial_condition=initial_condition, t=t_)
    sol2 = cloud_microphysics.FiniteCondensationTimeMicrophysics().integrate(initial_condition=initial_condition, t=t_)
