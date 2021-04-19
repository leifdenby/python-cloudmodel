from pyclouds.reference.atmos import stratification_profiles
from pyclouds.models import parcel as parcel_models
from pyclouds.models import microphysics as cloud_microphysics
from pyclouds import Var

import numpy as np

CloudModel = parcel_models.DryAirOnly


def test__integration_stop():
    """
    Check out that stopping criteria are working as expected
    """
    environment = stratification_profiles.Soong1973Dry()

    z_points = np.linspace(100.0, 10e3, 500)

    w0 = 0.1
    T0 = environment.temp(pos=0.0) + 0.2
    p0 = environment.p(pos=0.0)

    microphysics = cloud_microphysics.DummyMicrophysics()

    initial_condition = Var.make_state(
        r=500.0, w=w0, T=T0, q_v=0.0, q_l=0.0, q_r=0.0, q_i=0.0, p=p0
    )

    cloud_model = CloudModel(
        environment=environment, gamma=0.0, C_D=0.0, beta=0.2, microphysics=microphysics
    )

    cloud_profile = cloud_model.integrate(initial_condition, z_points)

    assert cloud_profile.z.max() > 0
