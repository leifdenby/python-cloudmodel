"""
Collection of microphysics routines for use with cloud-model integration.

return signature is:
    [dq_v, dq_l, dq_i,]

"""
import numpy as np

try:
    from ccfm.ccfmfortran import microphysics as ccfm_microphysics
    import testppp
except ImportError:
    # import pure python version instead
    from ccfm.ccfmpython import microphysics as ccfm_microphysics

def moist_adjustment(T, p, qv):
    """
    Adjust temperature and specific concentration of water vapour to saturation
    point, assuming pressure constants

    """
    T_new, qv_new = ccfm_microphysics.moist_adjust(tem=T, prs=p, q_v=qv)

    dq_v = qv_new - qv
    dq_l = qv - qv_new

    return [dq_v, dq_l, 0.0]
