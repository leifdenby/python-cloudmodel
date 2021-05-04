from pyclouds.reference.atmos import stratification_profiles, saturation_calculation
from pyclouds.models import parcel as parcel_models
from pyclouds.models import microphysics as microphysics_models
from pyclouds import Var
import pyclouds

import numpy as np


def test__profile_plots():
    """
    Integrate a cloud-model twice and plot the profiles
    """
    environment = stratification_profiles.Soong1973Dry()

    z_points = np.linspace(100.0, 4e3, 500)

    w0 = 0.1
    T0 = environment.temp(0.0) + 0.2
    p0 = environment.p(0.0)

    microphysics = microphysics_models.DummyMicrophysics()

    initial_condition = Var.make_state(
        r=500.0, w=w0, T=T0, q_v=0.0, q_l=0.0, q_r=0.0, q_i=0.0, p=p0
    )

    cloud_model = parcel_models.DryAirOnly(
        environment=environment, gamma=0.0, C_D=0.0, beta=0.0, microphysics=microphysics
    )
    profile_no_drag = cloud_model.integrate(
        initial_condition=initial_condition, z=z_points
    )

    cloud_model = parcel_models.DryAirOnly(
        environment=environment, gamma=0.0, C_D=1.0, beta=0.0, microphysics=microphysics
    )
    profile_with_drag = cloud_model.integrate(
        initial_condition=initial_condition, z=z_points
    )

    pyclouds.plot.parcel.plot_profiles([profile_no_drag, profile_with_drag])
