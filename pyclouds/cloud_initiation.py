"""
Contains a number of approaches for calculating the state a cloud-base which sets
the initial condition for cloud-profile integration.
"""

from common import default_constants, Var
from cloud_microphysics import MoistAdjustmentMicrophysics
import parameterisations

def compute_LCL(environment, F0):
    """
    Compute the height at which an adiabatically lifted sub-saturated parcel
    will condense given an environmental profile and ground-base conditions for the parcel.
    """
    # load constants
    cp_d = default_constants.get('cp_d')
    cp_v = default_constants.get('cp_v')
    R_d = default_constants.get('R_d')
    R_v = default_constants.get('R_v')
    g = default_constants.get('g')

    qv = F0[Var.q_v]
    ql = F0[Var.q_l]
    qr = F0[Var.q_r]
    qi = F0[Var.q_i]
    assert ql == 0.0 and qr == 0.0 and qi == 0.0
    qd = 1.0 - qv

    cp_g = qd*cp_d + qv*cp_v
    R_g = qd*R_d + qv*R_v
    T = F0[Var.T]

    # XXX: fixed integration step for now
    dz = 1.0
    z = 0.0
    while True:
        p = environment.p(z)
        qv_sat = parameterisations.pv_sat.qv_sat(T=T, p=p)

        Sw = qv/qv_sat

        if Sw > 1.0:
            break
        if z > 10e3:
            raise Exception()

        # temperature decreases at at dry adiabatic lapse rate
        dTdz = -g/cp_g

        z += dz
        T += dz*dTdz

    return z, T
