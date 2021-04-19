"""
Tests for integration of the cloud-equations will full thermodynamics
representation
"""
from pyclouds.reference.atmos import stratification_profiles
from pyclouds.models import parcel as parcel_models
from pyclouds.models import microphysics as cloud_microphysics
from pyclouds import Var

import numpy as np


CloudModel = parcel_models.FullThermodynamicsCloudEquations


def test__acceleration_by_latent_heat():
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
    # p, r, w, T, q_v, q_r, q_l, q_i

    beta = 0.0
    microphysics = cloud_microphysics.FiniteCondensationTimeMicrophysics()

    cloud_model = CloudModel(
        environment=environment,
        gamma=0.0,
        C_D=0.0,
        beta=beta,
        microphysics=microphysics,
    )

    initial_condition = Var.make_state(
        r=500.0, w=w0, T=T0, q_v=0.0, q_l=0.0, q_r=0.0, q_i=0.0, p=p0
    )
    profile_dry_parcel = cloud_model.integrate(
        initial_condition, z_points, fail_on_nan=True
    )

    initial_condition = Var.make_state(
        r=500.0, w=w0, T=T0, q_v=0.012, q_l=0.0, q_r=0.0, q_i=0.0, p=p0
    )
    profile_moist_parcel = cloud_model.integrate(
        initial_condition, z_points, fail_on_nan=True
    )

    assert profile_dry_parcel.z[-1] < profile_moist_parcel.z[-1]


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
    # p, r, w, T, q_v, q_r, q_l, q_i

    microphysics = cloud_microphysics.DummyMicrophysics()

    initial_condition = Var.make_state(
        r=500.0, w=w0, T=T0, q_v=0.0, q_l=0.0, q_r=0.0, q_i=0.0, p=p0
    )

    cloud_model = CloudModel(
        environment=environment, gamma=0.0, C_D=0.0, beta=0.0, microphysics=microphysics
    )
    profile_no_drag = cloud_model.integrate(initial_condition, z_points)

    cloud_model = CloudModel(
        environment=environment, gamma=0.0, C_D=1.0, beta=0.0, microphysics=microphysics
    )
    profile_with_drag = cloud_model.integrate(initial_condition, z_points)

    assert profile_with_drag.z[-1] < profile_no_drag.z[-1]


def test__entrainment():
    """
    With only dry air entrainment should still cause the parcel to slow down as colder air is being mixed in.
    """
    environment = stratification_profiles.Soong1973Dry()

    z_points = np.linspace(100.0, 4e3, 500)

    w0 = 0.1
    T0 = environment.temp(0.0) + 0.2
    p0 = environment.p(0.0)
    # p, r, w, T, q_v, q_r, q_l, q_i

    microphysics = cloud_microphysics.DummyMicrophysics()

    initial_condition = Var.make_state(
        r=500.0, w=w0, T=T0, q_v=0.0, q_l=0.0, q_r=0.0, q_i=0.0, p=p0
    )

    cloud_model_no_entrainment = CloudModel(
        environment=environment, gamma=0.0, C_D=0.0, beta=0.0, microphysics=microphysics
    )
    cloud_model_with_entrainment = CloudModel(
        environment=environment, gamma=0.0, C_D=0.0, beta=0.2, microphysics=microphysics
    )
    cloud_model_dry = parcel_models.DryAirOnly(
        environment=environment, gamma=0.0, D=0.0, beta=0.2, microphysics=microphysics
    )

    dFdz1 = cloud_model_no_entrainment.dFdz(initial_condition, 0.0)
    dFdz2 = cloud_model_with_entrainment.dFdz(initial_condition, 0.0)
    dFdz_dry = cloud_model_with_entrainment.dFdz(initial_condition, 0.0)

    # although only dry air is being entrained the effect of entrainment will
    # change the temperature as the entrained air will be cooler
    assert dFdz1[Var.T] > dFdz2[Var.T]

    # effect on momentum should be the same for dry and moist models
    assert abs(dFdz1[Var.w] - dFdz_dry[Var.w]) < 1.0e-4

    profile_with_entrainment__dry_model = cloud_model_dry.integrate(
        initial_condition, z_points
    )
    profile_with_entrainment = cloud_model_with_entrainment.integrate(
        initial_condition, z_points
    )
    profile_no_entrainment = cloud_model_no_entrainment.integrate(
        initial_condition, z_points
    )

    # with only dry air the max height should be near the same as with the dry
    # model. We relax the condition a bit since the buoyancy is formulated
    # differently in the two models
    dz = abs(profile_with_entrainment__dry_model.z[-1] - profile_with_entrainment.z[-1])
    assert dz / profile_with_entrainment__dry_model.z[-1] < 0.1

    # entraining ambient are with a dry plume should cause the plume to reach a lower height
    assert profile_with_entrainment.z[-1] < profile_no_entrainment.z[-1]


def test__finite_condensation_time_microphysics():
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
    # p, r, w, T, q_v, q_r, q_l, q_i

    beta = 0.2
    microphysics = cloud_microphysics.FiniteCondensationTimeMicrophysics()

    cloud_model = CloudModel(
        environment=environment,
        gamma=0.0,
        C_D=0.0,
        beta=beta,
        microphysics=microphysics,
    )

    initial_condition = Var.make_state(
        r=500.0, w=w0, T=T0, q_v=0.0, q_l=0.0, q_r=0.0, q_i=0.0, p=p0
    )

    # just run the integration for now so we can test that it doesn't fail
    cloud_model.integrate(initial_condition, z_points)
