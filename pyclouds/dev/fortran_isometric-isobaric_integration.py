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
T0 = 288.0
p0 = 88676.0  # [Pa]

initial_condition = Var.make_state(T=T0, p=p0, q_v=q_v, q_l=q_l)


t_ = np.linspace(0.0, 17.0, 100)

SolverClass = odespy.Euler

solutions = []
# solutions.append(cloud_microphysics.FiniteCondensationTimeMicrophysics(r_crit=6e-6, constants=um_constants).integrate(initial_condition=initial_condition, t=t_, p0=p0, SolverClass=SolverClass))
# solutions.append(cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants).integrate(initial_condition=initial_condition, t=t_, p0=p0, SolverClass=SolverClass))
# solutions.append(cloud_microphysics.FC_min_radius(constants=um_constants).integrate(initial_condition=initial_condition, t=t_, p0=p0, SolverClass=SolverClass))

model_constraint = "isobaric"
solutions.append(
    cloud_microphysics.ExplicitFortranModel(
        model_constraint=model_constraint
    ).integrate(initial_condition=initial_condition, t=t_, SolverClass=SolverClass)
)

model_constraint = "isometric"
solutions.append(
    cloud_microphysics.ExplicitFortranModel(
        model_constraint=model_constraint
    ).integrate(initial_condition=initial_condition, t=t_, SolverClass=SolverClass)
)

plot.ioff()
plot_hydrometeor_evolution(
    solutions,
    variables=[
        "q_v",
        "p",
        "q_l",
        "T",
    ],
    legend_loc="upper right",
)
plot.show()
plot.draw()
