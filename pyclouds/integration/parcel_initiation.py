"""
Contains a number of approaches for calculating the state a cloud-base which sets
the initial condition for cloud-profile integration.
"""

from .. import Var
from ..reference.constants import default_constants
from ..models.microphysics import MoistAdjustmentMicrophysics
from ..reference import parameterisations


class CloudbaseNotFoundException(Exception):
    pass


def compute_LCL(environment, F0):
    """
    Compute the height at which an adiabatically lifted sub-saturated parcel
    will condense given an environmental profile and ground-base conditions for the parcel.
    """
    # load constants
    cp_d = default_constants.get("cp_d")
    cp_v = default_constants.get("cp_v")
    R_d = default_constants.get("R_d")
    R_v = default_constants.get("R_v")
    g = default_constants.get("g")

    qv = F0[Var.q_v]
    ql = F0[Var.q_l]
    qr = F0[Var.q_r]
    qi = F0[Var.q_i]
    assert ql == 0.0 and qr == 0.0 and qi == 0.0
    qd = 1.0 - qv

    cp_g = qd * cp_d + qv * cp_v
    R_g = qd * R_d + qv * R_v
    T = F0[Var.T]

    # XXX: fixed integration step for now
    dz = 1.0
    z = 0.0
    while True:
        p = environment.p(z)
        qv_sat = parameterisations.pv_sat.qv_sat(T=T, p=p)

        Sw = qv / qv_sat

        if Sw > 1.0:
            return z, T
        if z > 10e3:
            raise CloudbaseNotFoundException

        # temperature decreases at at dry adiabatic lapse rate
        dTdz = -g / cp_g

        z += dz
        T += dz * dTdz

    raise CloudbaseNotFoundException


def original_CCFM_cloudbase(environment, dqv=0.0, dT=0.0):
    """
    Based of `cloudbase` in `mo_ccfm_cloudbase.f90`
    https://github.com/leifdenby/ccfm/blob/master/src/mo_ccfm_cloudbase.f90

    Returns a tuple of cloudbase (height, temperature, water vapour content)

      cl_base = 1
      cond    = .false.

      T_c = T_e(level)
      q_c = q_e(level)

      base: do l = level-1, 2 ,-1

         ! dry adiabatic ascent
         T_c =  ( cpd * T_c + geo(l+1) - geo(l) ) / cpd
         q_old = q_c

         ! adjustmend for moist ascent
         call moist_adjust( T_c, q_c, prs(l) )

         ! if condensation occurred
         if ( q_c < q_old ) cond = .true.

         buo =   T_c    * ( 1._dp + epsi * q_c )     &
               - T_e(l) * ( 1._dp + epsi * q_e(l) )

         if ( cond .and. (buo >= - basebuo) ) then
             cl_base = l
             exit
         end if

         if ( cond .and. (buo <= -(basebuo + 0.2_dp)) ) then
            exit
         endif

      end do base
    """
    mphys = MoistAdjustmentMicrophysics()

    basebuo = 0.5  # taken from `cloudbase.f90`

    cp_d = default_constants.get("cp_d")
    R_d = default_constants.get("R_d")
    R_v = default_constants.get("R_v")
    g = default_constants.get("g")

    epsi = R_d / R_v

    condensed = False

    dz = 10.0
    z = 0.0
    T_c = environment.temp(z) + dT
    qv_c = environment.q_v(z) + dqv

    # we don't expect a cloud-base above 4km
    while z < 4e3:
        z += dz
        T_c += -g / cp_d * dz
        p_c = environment.p(z)

        T_e = environment.temp(z)
        qv_e = environment.q_v(z)

        F = Var.make_state(T=T_c, p=p_c, q_v=qv_c)
        F_new = mphys._calc_adjusted_state(F, iterations=5)

        if F_new[Var.q_v] < F[Var.q_v]:
            condensed = True
            T_c = F_new[Var.T]

        buo = T_c * (1.0 + epsi * qv_c) - T_e * (1.0 + epsi * qv_e)

        if condensed and buo >= -basebuo:
            return z, T_c, qv_c

        if condensed and buo <= -(basebuo + 0.2):
            break

    raise CloudbaseNotFoundException
