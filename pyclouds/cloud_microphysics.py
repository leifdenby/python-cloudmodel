"""
Collection of microphysics routines for use with cloud-model integration.
"""
import numpy as np

try:
    from ccfm.ccfmfortran import microphysics as ccfm_microphysics
except ImportError:
    # import pure python version instead
    from ccfm.ccfmpython import microphysics as ccfm_microphysics


class MoistAdjustmentMicrophysics(object):
    """
    Adjust temperature and specific concentration of water vapour to saturation
    point, assuming pressure constants

    Uses `moist_adjust` method from `mo_ccfm_cloudbase.f90`
    """


    def __call__(self, F):
        # the fortran routine modifies variables in place, we don't want that
        T_new = np.array(T)
        qv_new = np.array(qv)

        mo_ccfm_cloudbase.moist_adjust(tem=T_new, prs=p, q_v=qv_new)

        dq_v = qv_new - qv

        Fs = np.copy(F.shape)

        Fs[Var.q_v] = qv_new
        Fs[Var.q_l] = F[Var.q_l] - dq_v

        return Fs
