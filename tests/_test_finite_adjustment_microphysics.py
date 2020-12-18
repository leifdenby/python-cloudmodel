import numpy as np

from pyclouds.models import microphysics as microphysics_models
from pyclouds import Var


def test_integration():
    initial_condition = Var.make_state(T=300.0, p=101325.0, q_v=0.0)

    t_ = np.linspace(0.0, 10.0, 100)

    microphysics_model = microphysics_models.MoistAdjustmentMicrophysics()
    sol = microphysics_model.integrate(initial_condition=initial_condition, t=t_)
    sol2 = microphysics_models.FiniteCondensationTimeMicrophysics().integrate(
        initial_condition=initial_condition, t=t_
    )
