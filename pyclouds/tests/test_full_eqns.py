from pycfd.reference.atmospheric_flow import stratification_profiles, saturation_calculation
from pyclouds import cloud_equations, plotting, cloud_microphysics
from pyclouds.common import Var

import numpy as np


CloudModel = cloud_equations.FullThermodynamicsCloudEquations

def test__acceleration_by_latent_heat():
    """
    Without entrainment of ambient air a moist parcel will rise much further
    than a dry parcel if the moist parcel reaches saturation, as the moist
    parcel will warm up from latent heat release.
    """
    environment = stratification_profiles.Soong1973Dry()

    z_points = np.linspace(100., 4e3, 500)

    w0 = 0.1
    T0 = environment.temp(0.0) + 0.2
    # p, r, w, T, q_v, q_r, q_l, q_i

    beta = 0.0
    microphysics = cloud_microphysics.MoistAdjustmentMicrophysics()

    cloud_model = CloudModel(environment=environment, gamma=0.0, D=0.0, beta=beta, microphysics=microphysics)
    
    initial_condition = [500.0, w0, T0, 0.0, 0., 0., 0.]
    profile_dry_parcel = cloud_model.integrate(initial_condition, z_points)

    initial_condition = [500.0, w0, T0, 0.012, 0., 0., 0.]
    profile_moist_parcel = cloud_model.integrate(initial_condition, z_points)

    assert profile_dry_parcel.z[-1] < profile_moist_parcel.z[-1]

def test__drag():
    """
    Without entrainment of ambient air a moist parcel will rise much further
    than a dry parcel if the moist parcel reaches saturation, as the moist
    parcel will warm up from latent heat release.
    """
    environment = stratification_profiles.Soong1973Dry()

    z_points = np.linspace(100., 4e3, 500)

    w0 = 0.1
    T0 = environment.temp(0.0) + 0.2
    # p, r, w, T, q_v, q_r, q_l, q_i

    microphysics = cloud_microphysics.DummyMicrophysics()
    
    initial_condition = [500.0, w0, T0, 0.0, 0., 0., 0.]

    cloud_model = CloudModel(environment=environment, gamma=0.0, D=0.0, beta=0.0, microphysics=microphysics)
    profile_no_drag = cloud_model.integrate(initial_condition, z_points)

    cloud_model = CloudModel(environment=environment, gamma=0.0, D=1.0, beta=0.0, microphysics=microphysics)
    profile_with_drag = cloud_model.integrate(initial_condition, z_points)

    assert profile_with_drag.z[-1] < profile_no_drag.z[-1]

def test__entrainment():
    """
    With only dry air entrainment should still cause the parcel to slow down as colder air is being mixed in.
    """
    environment = stratification_profiles.Soong1973Dry()

    z_points = np.linspace(100., 4e3, 500)

    w0 = 0.1
    T0 = environment.temp(0.0) + 0.2
    # p, r, w, T, q_v, q_r, q_l, q_i

    microphysics = cloud_microphysics.DummyMicrophysics()
    
    initial_condition = [500.0, w0, T0, 0.0, 0., 0., 0.]


    cloud_model_with_entrainment = CloudModel(environment=environment, gamma=0.0, D=0.0, beta=0.2, microphysics=microphysics)
    cloud_model_no_entrainment = CloudModel(environment=environment, gamma=0.0, D=0.0, beta=0.0, microphysics=microphysics)
    cloud_model_dry = cloud_equations.DryAirOnly(environment=environment, gamma=0.0, D=0.0, beta=0.2, microphysics=microphysics)

    dFdz1 = cloud_model_no_entrainment.dFdz(initial_condition, 0.0)
    dFdz2 = cloud_model_with_entrainment.dFdz(initial_condition, 0.0)
    dFdz_dry = cloud_model_with_entrainment.dFdz(initial_condition, 0.0)

    # since only dry air is being entrained the effect of entrainment should
    # not effect the temperature, since the moist static energy of the
    # environment and the cloud will be the same
    assert dFdz1[Var.T] == dFdz2[Var.T]

    # effect on moment should be the same for dry and moist models
    assert abs(dFdz1[Var.w] - dFdz_dry[Var.w]) < 1.e-4

    profile_with_entrainment__dry_model = cloud_model_dry.integrate(initial_condition, z_points)
    profile_with_entrainment = cloud_model_with_entrainment.integrate(initial_condition, z_points)
    profile_no_entrainment = cloud_model_no_entrainment.integrate(initial_condition, z_points)

    # with only dry air the max height should be near the same as with the dry
    # model. We relax the condition a bit since the buoyancy is formulated
    # differently in the two models
    dz = abs(profile_with_entrainment__dry_model.z[-1] - profile_with_entrainment.z[-1])
    assert  dz/profile_with_entrainment__dry_model.z[-1] < 0.1

    # entraining ambient are with a dry plume should cause the plume to reach a lower height
    assert profile_with_entrainment.z[-1] < profile_no_entrainment.z[-1]
