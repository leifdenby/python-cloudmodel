from pyclouds.cloud_initiation import compute_LCL
from pyclouds.common import Var
from pyclouds import parameterisations

from pycfd.reference.atmospheric_flow import stratification_profiles


def test_LCL():
    environment = stratification_profiles.Soong1973()

    dT = 0.1

    T0 = environment.temp(0.0) + dT
    p0 = environment.p(0.0)
    qv_sat = parameterisations.pv_sat.qv_sat(T=T0, p=p0)
    q_v = 0.6*qv_sat

    F0 = Var.make_state(T=T0, p=p0, q_v=q_v)

    z_clb, T = compute_LCL(environment=environment, F0=F0)

    assert z_clb > 0.0
    assert z_clb < 1300.


    # higher initial temperature should increase cloud-base height
    dT = 1.0
    T0 = environment.temp(0.0) + dT
    p0 = environment.p(0.0)
    qv_sat = parameterisations.pv_sat.qv_sat(T=T0, p=p0)
    q_v = 0.6*qv_sat

    F0 = Var.make_state(T=T0, p=p0, q_v=q_v)

    z_clb2, T = compute_LCL(environment=environment, F0=F0)


    assert z_clb < z_clb2


    # higher initial RH should lower cloud-base height
    dT = 0.1
    T0 = environment.temp(0.0) + dT
    p0 = environment.p(0.0)
    qv_sat = parameterisations.pv_sat.qv_sat(T=T0, p=p0)
    q_v = 0.8*qv_sat

    F0 = Var.make_state(T=T0, p=p0, q_v=q_v)

    z_clb3, T = compute_LCL(environment=environment, F0=F0)


    assert z_clb > z_clb3
