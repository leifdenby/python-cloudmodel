"""
Utility script for creating vertical profiles of convective clouds using input
profiles (`INPUT_profile`) from ATHAM and a particular cloud-profile
"""
import argparse
import numpy as np
import matplotlib.pyplot as plot

import pyclouds.cloud_equations
from pyclouds.common import Var
from pyclouds import cloud_initiation, parameterisations, cloud_microphysics
import pyclouds.plotting

from pycfd.reference.atmospheric_flow import stratification_profiles

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("model")
parser.add_argument("-N", help="number of clouds", default=5, type=int)
parser.add_argument("-w_c", help="cloud-base velocity", default=0.1)
parser.add_argument("-r_c__max", help="cloud-base max_radius", default=None, type=float)
parser.add_argument("-z_max", help="max integration height", default=10e3)
parser.add_argument(
    "-dT_b", help="bubble temperature perturbation", default=0.1, type=float
)
parser.add_argument("-dqv_b", help="bubble moisture perturbation", default=0.0)
parser.add_argument("-beta", help="entrainment coefficient", default=0.1, type=float)
parser.add_argument("-debug", help="debug", default=True)

args = parser.parse_args()
debug = args.debug

# environment
# environment = stratification_profiles.Soong1973()
environment = stratification_profiles.RICO()

# calculate cloud-base height and temperature for LCL
dT_b = args.dT_b
T0 = environment.temp(0.0) + dT_b
p0 = environment.p(0.0)

qv_sat0 = parameterisations.pv_sat.qv_sat(T=T0, p=p0)
dqv_b = args.dqv_b
qv0 = environment.rel_humidity(0.0) * qv_sat0 + dqv_b

F0_bubble = Var.make_state(T=T0, p=p0, q_v=qv0)
z_clb, T_clb = cloud_initiation.compute_LCL(environment=environment, F0=F0_bubble)

if debug:
    print(
        """Integrating with
    z_clb, T_clb = {z_clb}m, {T_clb}T
    """.format(
            z_clb=z_clb, T_clb=T_clb
        )
    )

# setup cloud model
model_name = args.model
ModelClass = getattr(pyclouds.cloud_equations, model_name)
microphysics = cloud_microphysics.FiniteCondensationTimeMicrophysics()

beta = args.beta
model = ModelClass(
    environment=environment, gamma=0.2, beta=beta, D=0.1, microphysics=microphysics
)
z_max = args.z_max
N_clouds = args.N
z = np.linspace(z_clb, args.z_max, (z_max - z_clb) / 10.0)

# max cloud radius
r_max = args.r_c__max
if r_max is None:
    r_max = z_clb

integrator = pyclouds.cloud_equations.odespy.RKFehlberg
profiles = []
for r_c in np.linspace(100.0, r_max, N_clouds, endpoint=True):
    # make cloud-base state
    w_c = args.w_c
    p0 = environment.p(z_clb)
    F0_cldbase = Var.make_state(T=T_clb, p=p0, r=r_c, w=w_c, q_v=qv0)

    # integrate
    profile = model.integrate(initial_condition=F0_cldbase, z=z, SolverClass=integrator)
    profiles.append(profile)

plot.ioff()
label_f = lambda profile: "$r_c={r_c}m$".format(r_c=profile.F[0, Var.r])
pyclouds.plotting.plot_profiles(
    profiles, ["r", "w", "T", "q_v", "q_l", "q_r"], label_f=label_f, labels_ncol=6
)
title = """Integration from LCL in {environment}
model: {model}
$dT_{{bubble}}={dT_b}K$, $dq_{{v,bubble}}={dqv_b}kg/kg$, $z_{{clb}}={z_clb}m$
integrator: {integrator}
""".format(
    environment=environment,
    model=str(model),
    dT_b=dT_b,
    dqv_b=dqv_b,
    integrator=str(integrator),
    z_clb=z_clb,
)
plot.suptitle(title)
plot.show()
