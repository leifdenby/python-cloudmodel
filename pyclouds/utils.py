import numpy as np

from common import Var

default_constants = {
    # Heat capacities
    "cp_d": 1004.64, # J/kg/K
    "cv_d": 717.60, # J/kg/K
    "cp_v": 1874.0, # J/kg/K
    "cv_v": 1402.5, # J/kg/K
    # Teten's formula for saturation vapour pressure (constants from Bechtold's notes)
    # over liquid water
    "p0vs": 611.2,  # [Pa]
    "a0_lq": 17.67,
    "a1_lq": -32.19,
    "a0_ice": 22.587,
    "a1_ice": 0.7,
}


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

class Utils:
    def __init__(self, constants):
        self.constants = AttrDict(constants)
        if not hasattr(constants, 'R_d'):
            self.constants.R_d = self.constants.cp_d - self.constants.cv_d
            self.constants.R_v = self.constants.cp_v - self.constants.cv_v

    def pv_sat_liquid(self, T):
        p0vs = self.constants.p0vs
        a0_lq = self.constants.a0_lq
        a1_lq = self.constants.a1_lq

        return p0vs*np.exp((a0_lq*(T-273.15)/(T+a1_lq)))

    def pv_sat_ice(self, T):
        p0vs = self.constants.p0vs
        a0_ice = self.constants.a0_ice
        a1_ice = self.constants.a1_ice

        return p0vs*np.exp((a0_ice*(T-273.15)/(T+a1_ice)))

    def pv_sat(self, T):
        v = np.zeros(np.array(T).shape)

        if v.shape == ():
            if T > 273.15:
                return self.pv_sat_liquid(T)
            else:
                return self.pv_sat_ice(T)
        else:
            T = np.array(T)
            idx_liquid = np.array(T) > 273.15
            idx_ice = idx_liquid == False
            v[idx_liquid] = self.pv_sat_liquid(T[idx_liquid])
            v[idx_ice] = self.pv_sat_ice(T[idx_ice])
            return v

    def qv(self, T, p, pv):
        """
        Calculate the specific concentration of water vapour given the partial
        pressure on the water phase.
        """
        # TODO: Does this equation assume only water vapour and dry air is present?

        R_d = self.constants.R_d
        R_v = self.constants.R_v

        epsilon = R_d/R_v

        return (epsilon*pv)/(p-(1-epsilon)*pv)

    def qv_sat(self, T, p):
        pv = self.pv_sat(T)
        return self.qv(T=T, p=p, pv=pv)

# pv_sat = lambda T: 611.2*np.exp((17.67*(T-273.15)/(T-32.19)))


    def moist_static_energy(self, F, z):
        q_v = F[...,Var.q_v]
        q_l = F[...,Var.q_l]
        q_i = F[...,Var.q_i]
        T = F[...,Var.T]

        g = self.constants.g
        L_v = self.constants.L_v
        cp_d = self.constants.cp_d
        cp_v = self.constants.cp_v
        L_s = self.constants.L_s
        L_v = self.constants.L_s

        q_d = 1.0 - q_v - q_l - q_i

        cp_m = cp_d*q_d + cp_v*(q_v + q_l + q_i)

        return cp_m*T - q_l*L_v - q_i*L_s + g*z

_utils_instance = None
def qv_sat(T, p):
    global _utils_instance
    if _utils_instance is None:
        print "Using default constants"
        _utils_instance = Utils(constants=default_constants)

    return _utils_instance.qv_sat(T, p)



# dry air and water vapour EoS:
# cp = lambda qv: (1-qv)*cp_d + qv*cp_v
# cv = lambda qv: (1-qv)*cv_d + qv*cv_v
# R = lambda qv: cp(qv) - cv(qv)
# kappa = lambda qv: R(qv)/cp(qv)

