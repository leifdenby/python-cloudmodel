import odespy
import numpy as np

default_constants = {
    "R_d": 287.05,
    "R_v": 461.51,
    "L_v": 2.5008e6,
    "L_s": 2.8345e6,
    "cp_d": 1005.46,
    "cv_d": 717.60, # J/kg/K
    "cp_v": 1859.0,
    "cv_v": 1402.5, # J/kg/K
    "cp_l": 4183.0,
    "cv_l": 4183.0, # same as cp as liquid is assumed incompressible
    "rho_l": 1000.,
    "rho_i": 500.,
    "g": 9.80665,
    "pv_sat": {  # constants for calculating saturation vapour pressure with Teten's formula
        "p0vs": 611.2,
        "a0_lq": 17.67,
        "a1_lq": -32.19,
    }
}

CCFM_constants = {
    "R_d": 287.05,  # 'Rd'
    "R_v": 461.51,  # 'Rv'
    "L_v": 2.5008e6,  # 'aLv'
    "L_s": 2.8345e6,  # 'aLs'
    "cp_d": 1005.46,  # 'Cpd'
    "cp_v": 1859.0,  # not specified in CCFM?
    "rho_l": 1000.,  # 'Rho_w'
    "rho_i": 500.,  # 'Rho_i'
    "g": 9.80665,  # 'g'
    "pv_sat": {
        'p0vs': 610.78,
        'a0_lq': 17.269,
        'a1_lq': -35.86,
        'a0_ice': 21.875,
        'a1_ice': -7.66,
    }
}

ATHAM_constants = {
    "cp_d": 1004.64,
    "cv_d": 717.60,
    "cp_v": 1870.,
    "cv_v": 1402.5,
    "L_v": 2500800.,
    "cp_l": 4183.,
    "cp_i": 2103.,
    "rho_l": 1000.,
    "rho_i": 917.,
    "pv_sat": {
        'p0vs': 610.7,
        'a0_lq': 17.25,
        'a1_lq': -36.0,
        'a0_ice': 22.33,
        'a1_ice': -2.,
    },
}

def make_related_constants(constants):
    if 'cp_d' in constants and 'cp_v' in constants and not 'R_d' in constants:
        constants['R_d'] = constants['cp_d'] - constants['cv_d']
        constants['R_v'] = constants['cp_v'] - constants['cv_v']
    return constants

um_constants = {
    "cp_d": 1004.64,
    "cv_d": 717.60,
    "cp_v": 1864.,
    "cv_v": 1402.55,
    "cp_l": 4183.,
    "cv_l": 4183.,
    "cp_i": 2103.,
    "cv_i": 2103.,
    "L_v": 2500.8,
    "rho_l": 1000.,
    "pv_sat": {  # constants for calculating saturation vapour pressure with Teten's formula
        "p0vs": 610.7,
        "a0_lq": 17.25,
        "a1_lq": 36.
    },
}

# r, w, T, q_v, q_r, q_l, q_i
class Var:
    r = 0
    w = 1
    T = 2
    q_v = 3
    q_r = 4
    q_l = 5
    q_i = 6
    z = 7
    p = 8

    names = ['r', 'w', 'T', 'q_v', 'q_r', 'q_l', 'q_i', 'z', 'p',]
    NUM = len(names)

    @staticmethod
    def print_formatted(v, formatting='%g'):
        print ",\t".join([("%s=" + formatting) % (Var.names[i], v[i]) for i in range(Var.NUM)])

    @staticmethod
    def repr(v, formatting='%g'):
        units = { 'w': 'm/s', 'r': 'm', 'T': 'K', 'z': 'm', 'p': 'Pa',}
        return ", ".join([r"$%s=%g%s$" % (Var.names[i], v[i], units.get(Var.names[i], '')) for i in range(Var.NUM)])

    @staticmethod
    def make_state(**kwargs):
        s = np.zeros((Var.NUM))
        for k, v in kwargs.items():
            s[getattr(Var, k)] = v

        return s

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self



REQUIRED_POSITIVE = np.zeros((Var.NUM))
# REQUIRED_POSITIVE[Var.r] = 1
# REQUIRED_POSITIVE[Var.q_v] = 1
REQUIRED_POSITIVE[Var.q_l] = 1
# REQUIRED_POSITIVE[Var.q_r] = 1
# REQUIRED_POSITIVE[Var.q_i] = 1

class MySolver(odespy.RKFehlberg):

    def advance(self):
        f, n, rtol, atol = self.f, self.n, self.rtol, self.atol
        u_n, t_n, t_np1 = self.u[n], self.t[n], self.t[n+1]
        dt = t_np1 - t_n

        dt_max__all = np.abs(f(u_n, t_n)/u_n)
        dt_max__all[REQUIRED_POSITIVE == 0] = 1.0e31

        self.max_step = min(np.min(dt_max__all), self.first_step)
        print self.max_step, dt_max__all
        self.first_step = self.max_step
        import ipdb
        with ipdb.launch_ipdb_on_exception():
            return odespy.RKFehlberg.advance(self)
