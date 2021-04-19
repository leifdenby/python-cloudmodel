# coding: utf-8
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plot
import numpy as np
import odespy
import sys

from pyclouds import cloud_microphysics, parameterisations
from pyclouds.common import Var, make_related_constants
from pyclouds.plotting import plot_hydrometeor_evolution

from unified_microphysics.tests.test_common import um_constants


constants = make_related_constants(um_constants)
t_max = 5.0

for dt in [t_max / 50, t_max / 100, t_max / 200]:
    t_ = np.linspace(0.0, 5.0, t_max / dt)
    for q_l in [1.0e-3, 1.0e-4, 1.0e-5]:

        initial_condition = np.zeros((Var.NUM))
        q_v = 0.0
        T = 288.0
        p0 = 88676.0  # [Pa]

        initial_condition[Var.q_v] = q_v
        initial_condition[Var.q_l] = q_l
        initial_condition[Var.T] = T

        solutions = []
        solutions.append(
            cloud_microphysics.OldATHAMKesslerFortran().integrate(
                initial_condition=initial_condition, t=t_, p0=p0
            )
        )
        solutions.append(
            cloud_microphysics.FiniteCondensationTimeMicrophysics(
                constants=constants
            ).integrate(initial_condition=initial_condition, t=t_, p0=p0)
        )

        plot_hydrometeor_evolution(
            solutions,
            variables=[
                "q_l",
                "T",
            ],
            legend_loc="upper right",
        )

        L_v = um_constants.get("L_v")
        cp_d = um_constants.get("cp_d")

        dT_approx = L_v / cp_d * q_l

        title = r"Initial condition: $q_v=0.0kg/kg$, $q_l=%gkg/kg$" % q_l
        title += (
            "\n"
            + r"with $c_{p,d}=%g\frac{J}{kg K}$, $L_v=%gJ/kg$ expect $\Delta T\approx \frac{L_v q_l}{c_{p,d}}=%gK$"
            % (cp_d, L_v, dT_approx)
        )
        title += "\n" + r"$\Delta t=%gs$" % dt

        plot.suptitle(title)
        plot.savefig("old-kessler_vs_finite-evapcond-rate_ql-%f_dt-%fs.png" % (q_l, dt))

# In[ ]:
