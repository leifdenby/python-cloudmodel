from .fortranlib import mo_ccfm_cloudbase

import numpy as np


def moist_adjust(tem, prs, q_v):
    """
    Wrapper for CCFM routine to compute moist adjustment at constant pressure.
    To avoid inplace update of state variables.

    Uses `moist_adjust` method from `mo_ccfm_cloudbase.f90`
    """

    # the fortran routine modifies variables in place, we don't want that
    T_new = np.array(tem)
    qv_new = np.array(q_v)

    mo_ccfm_cloudbase.moist_adjust(tem=T_new, prs=prs, q_v=qv_new)

    return T_new, qv_new
