from pyclouds.reference.atmos import stratification_profiles
from pyclouds.models import parcel as parcel_models
from pyclouds.models import microphysics as cloud_microphysics
from pyclouds import Var

import numpy as np


def test__integration_stop_altitude_too_high():
    """
    Check out that stopping criteria are working as expected
    """
    environment = stratification_profiles.Soong1973Dry()

    z_points = np.linspace(100.0, 100e3, 500)

    w0 = 0.1
    T0 = environment.temp(pos=0.0) + 0.2
    p0 = environment.p(pos=0.0)

    microphysics = cloud_microphysics.DummyMicrophysics()

    initial_condition = Var.make_state(
        r=500.0, w=w0, T=T0, q_v=0.0, q_l=0.0, q_r=0.0, q_i=0.0, p=p0
    )

    cloud_model = parcel_models.FixedRiseRateParcel(
        environment=environment, beta=0.2, microphysics=microphysics, w0=w0
    )

    cloud_profile = cloud_model.integrate(initial_condition, z_points)

    assert cloud_profile.z.max() > 0
    assert cloud_profile.integration_stopping_reason
    assert "height_unphysical" in cloud_profile.integration_stopping_reason
