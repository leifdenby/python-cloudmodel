"""
Collection of cloud-model equations.
"""
import odespy

# r, w, T, q_v, q_r, q_l, q_i
class Var:
    r = 0
    w = 1
    T = 2
    q_v = 3
    q_r = 4
    q_l = 5
    q_i = 6

class CloudModel(object):
    def __init__(self, environment, constants):
        self.T_e = environment.get('T_e')
        self.p_e = environment.get('p_e')
        self.g = constants.get('g')
        self.R_d = constants.get('R_d')
        self.cp_d = constants.get('cp_d')

    def dFdz(F, z):
        raise NotImplemented

    def integrate(self, initial_condition, z, SolverClass=odespy.RKFehlberg):
        solver = SolverClass(self.dFdz, rtol=0.0, atol=1e-6,)
        solver.set_initial_condition(initial_condition)

        stopping_criterion = lambda F, z, k: F[k,1] <= 0.0 or F[k,2] > 300.

        F, z = solver.solve(z, stopping_criterion)

        F = F[:-2]
        z = z[:-2]

        return F, z

class NoMicrophysicsNoEntrainment:
    mu = 0.0  # entrainment rate [1/m]
    gamma = 0.5  # virtual mass coefficient [1]
    D = 0.0  # drag coefficient [1]
    Lv = 2.5008e6  # latent heat of vapourisation [J/kg]

    def dwdz(z, r, w, T):
        """
        Also requires: environment temperature (density), so that buoyancy can be computed
        """

        B = (T-Te(z))/T

        return 1./w*(g/(1+gamma)*B - mu*w**2 - D*w**2/r)

    def dTdz(z, r, w, T):
        Te_ = Te(z)
        p_ = p(z)

        qsat_w = saturation_calculation.qv_sat(T=T, p=p_)

        return (-g/cp_d*(1+(Lv * qsat_w)/(R_d*T)) - mu*(Te_-T))

    def drdz(z, r, w, T, dwdz_, dTdz_):
        #dwdz_ = dwdz(z, r, w, T)
        #dTdz_ = dTdz(z, r, w, T)

        return 2./r *( (g/(R_d*T) + 1./T*dTdz_) - 1./w * dwdz_ + mu/r)

    def dFdz(F, z):
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]

        dwdz_ = dwdz(z, r, w, T)
        dTdz_ = dTdz(z, r, w, T)
        drdz_ = drdz(z, r, w, T, dwdz_, dTdz_)

        return [drdz_, dwdz_, dTdz_, 0., 0., 0., 0.,]

class Wagner2009(CloudModel):
    """
    Moist saturation cloud-model, from Till's PhD thesis
    """
    Lv = 2.5008e6  # latent heat of vapourisation [J/kg]

    def __init__(self, **kwargs):
        """
        mu: entrainment rate [1/m]
        gamma: virtual mass coefficient [1]
        D: drag coefficient [1]
        """
        self.mu = kwargs.pop('mu', 0.0)
        self.gamma = kwargs.pop('gamma', 0.5)
        self.D = kwargs.pop('D', 0.0)

        self.microphysics = kwargs.pop('microphysics')

        self.qv_e = kwargs['environment'].get('qv_e')
        super(Wagner2009, self).__init__(**kwargs)

    def dwdz(self, z, r, w, T):
        """
        Also requires: environment temperature (density), so that buoyancy can be computed
        """
        g, gamma, mu, D = self.g, self.gamma, self.mu, self.D

        B = (T-self.T_e(z))/T

        return 1./w*(g/(1+gamma)*B - mu*w**2 - D*w**2/r)

    def dTdz__dQdz(self, z, r, w, T, q_v):
        g, cp_d, mu = self.g, self.cp_d, self.mu

        Te_ = self.T_e(z)
        p_e = self.p_e(z)

        dTdz_s = -g/cp_d - mu*(T - Te_)

        # Intermediate values
        # TODO: entrainment of air
        T_s = T + dTdz_s
        qv__s = q_v

        dQ = self.microphysics(T=T_s, p=p_e, qv=qv__s)
        dq_v = dQ[0]

        dTdz__latent_heat = self.Lv*dq_v

        dTdz_ = dTdz_s + dTdz__latent_heat

        return dTdz_, dQ

    def drdz(self, z, r, w, T, dwdz_, dTdz_):
        g, R_d, mu = self.g, self.R_d, self.mu

        return r/2. *( (g/(R_d*T) + 1./T*dTdz_) - 1./w * dwdz_ + mu/r)

    def dFdz(self, F, z):
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]
        q_v = F[Var.q_v]

        dwdz_ = self.dwdz(z, r, w, T)
        dTdz_, dQdz_ = self.dTdz__dQdz(z, r, w, T, q_v)
        drdz_ = self.drdz(z, r, w, T, dwdz_, dTdz_)

        # r, w, T, q_v, q_r, q_l, q_i
        return [drdz_, dwdz_, dTdz_, dQdz_[0], dQdz_[1], 0., 0.,]

class DryAirOnly(CloudModel):
    def __init__(self, **kwargs):
        """
        mu: entrainment rate [1/m]
        gamma: virtual mass coefficient [1]
        D: drag coefficient [1]
        """
        self.mu = kwargs.pop('mu', 0.0)
        self.gamma = kwargs.pop('gamma', 0.5)
        self.D = kwargs.pop('D', 0.0)
        super(DryAirOnly, self).__init__(**kwargs)

    def dwdz(self, z, r, w, T):
        """
        Also requires: environment temperature (density), so that buoyancy can be computed
        """
        g, gamma, mu, D = self.g, self.gamma, self.mu, self.D

        B = (T-self.T_e(z))/T

        return 1./w*(g/(1+gamma)*B - mu*w**2 - D*w**2/r)

    def dTdz(self, z, r, w, T):
        g, cp_d, mu = self.g, self.cp_d, self.mu

        Te_ = self.T_e(z)

        return -g/cp_d - mu*(T - Te_)

    def drdz(self, z, r, w, T, dwdz_, dTdz_):
        g, R_d, mu = self.g, self.R_d, self.mu

        return r/2. *( (g/(R_d*T) + 1./T*dTdz_) - 1./w * dwdz_ + mu/r)

    def dFdz(self, F, z):
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]

        dwdz_ = self.dwdz(z, r, w, T)
        dTdz_ = self.dTdz(z, r, w, T)
        drdz_ = self.drdz(z, r, w, T, dwdz_, dTdz_)

        return [drdz_, dwdz_, dTdz_, 0., 0., 0., 0.,]
