
# coding: utf-8

# In[1]:

import matplotlib.pyplot as plot
import numpy as np
import odespy


from pyclouds import cloud_microphysics, parameterisations
from pyclouds.common import Var
from pyclouds.plotting import plot_hydrometeor_evolution

from unified_microphysics.tests.test_common import um_constants


# In[ ]:

q_v = 1.1e-2
q_l = 2.0e-4
T0 = 288.
p0 = 88676.  # [Pa]

initial_condition = Var.make_state(T=T0, p=p0, q_v=q_v, q_l=q_l)


t_ = np.linspace(0., 17., 100)

SolverClass = odespy.Euler

def test_isobaric_isometric():
    model_constraint = 'isobaric'
    sol1 = cloud_microphysics.ExplicitFortranModel(model_constraint=model_constraint).integrate(initial_condition=initial_condition, t=t_, SolverClass=SolverClass)

    model_constraint = 'isometric'
    sol2 = cloud_microphysics.ExplicitFortranModel(model_constraint=model_constraint).integrate(initial_condition=initial_condition, t=t_, SolverClass=SolverClass)

    # pressure should be unchanged for isobaric integration
    assert sol1.F[-1][Var.p] == p0

    Sw = initial_condition[Var.q_v]/parameterisations.pv_sat.qv_sat(T=T0, p=p0)

    if Sw > 1.0:
        assert sol2.F[-1][Var.q_l] > q_l
        assert sol2.F[-1][Var.p] > p0
    else:
        assert sol2.F[-1][Var.q_l] < q_l
        assert sol2.F[-1][Var.p] < p0
