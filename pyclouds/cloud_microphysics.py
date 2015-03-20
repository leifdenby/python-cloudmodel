"""
Collection of microphysics routines for use with cloud-model integration.

return signature is:
    [dq_v, dq_l, dq_i,]

"""
import numpy as np

from pyclouds.ccfm.cloudbase import mo_ccfm_cloudbase

def moist_adjustment(T, p, qv):
    """
    Adjust temperature and specific concentration of water vapour to saturation
    point, assuming pressure constants

    Uses `moist_adjust` method from `mo_ccfm_cloudbase.f90`
    """

    # the fortran routine modifies variables in place, we don't want that
    T_new = np.array(T)
    qv_new = np.array(qv)

    mo_ccfm_cloudbase.moist_adjust(tem=T_new, prs=p, q_v=qv_new)

    dq_v = qv_new - qv
    dq_l = qv - qv_new

    return [dq_v, dq_l, 0.0]
