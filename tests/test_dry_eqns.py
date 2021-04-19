from pyclouds.reference.atmos import stratification_profiles, saturation_calculation
from pyclouds.models import parcel as parcel_models
from pyclouds.models import microphysics as microphysics_models
from pyclouds import Var

import numpy as np


def test__drag():
    """
    Without entrainment of ambient air a moist parcel will rise much further
    than a dry parcel if the moist parcel reaches saturation, as the moist
    parcel will warm up from latent heat release.
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

    assert profile_with_drag.z[-1] < profile_no_drag.z[-1]


def test__entrainment():
    """
    With only dry air entrainment should still cause the parcel to slow down as colder air is being mixed in.
    """
    CloudModel = parcel_models.DryAirOnly
    environment = stratification_profiles.Soong1973Dry()

    z_points = np.linspace(100.0, 4e3, 500)

    w0 = 0.1
    T0 = environment.temp(0.0) + 0.2
    p0 = environment.p(0.0)
    # p, r, w, T, q_v, q_r, q_l, q_i

    microphysics = microphysics_models.DummyMicrophysics()

    initial_condition = Var.make_state(
        r=500.0, w=w0, T=T0, q_v=0.0, q_l=0.0, q_r=0.0, q_i=0.0, p=p0
    )

    cloud_model = CloudModel(
        environment=environment, gamma=0.0, C_D=0.0, beta=0.0, microphysics=microphysics
    )
    profile_no_entrainment = cloud_model.integrate(initial_condition, z_points)

    cloud_model = CloudModel(
        environment=environment, gamma=0.0, C_D=0.0, beta=0.2, microphysics=microphysics
    )
    profile_with_entrainment = cloud_model.integrate(initial_condition, z_points)

    assert profile_with_entrainment.z[-1] < profile_no_entrainment.z[-1]
