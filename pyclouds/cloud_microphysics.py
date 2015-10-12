"""
Collection of microphysics routines for use with cloud-model integration.
"""
import numpy as np

from pyclouds.cloud_equations import Var

try:
    from ccfm.ccfmfortran import microphysics as ccfm_microphysics
except ImportError:
    # import pure python version instead
    from ccfm.ccfmpython import microphysics as ccfm_microphysics


class DummyMicrophysics(object):
    """
    Dummy microphysics implementation that doesn't affect the state at all.
    """


    def __call__(self, F, p):
        return F

class MoistAdjustmentMicrophysics(object):
    """
    Adjust temperature and specific concentration of water vapour to saturation
    point, assuming pressure constants

    Uses `moist_adjust` method from `mo_ccfm_cloudbase.f90`
    """


    def __call__(self, F, p):
        T = F[Var.T]
        qv = F[Var.q_v]

        T_new, qv_new = ccfm_microphysics.moist_adjust(tem=T, prs=p, q_v=qv)

        dq_v = qv_new - qv

        Fs = np.copy(F)

        Fs[Var.q_v] = qv_new
        Fs[Var.q_l] = F[Var.q_l] - dq_v
        Fs[Var.T] = T_new

        return Fs
