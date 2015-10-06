"""
Collection of cloud-model equations.
"""
import odespy
import warnings

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

        stopping_criterion = lambda F, z, k: F[k,1] <= 0.0 or F[k,2] > 300. or F[k,Var.T] < 0.0

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
    beta = 0.2  # entraiment-rate coefficient [1]
    gamma = 0.5  # virtual-mass coefficient [1]

    def __init__(self, **kwargs):
        self.microphysics = kwargs.pop('microphysics')
        self.qv_e = kwargs['environment'].get('qv_e')
        super(Wagner2009, self).__init__(**kwargs)

    def mu(self, r):
        return self.beta/r

    def dwdz(self, z, r, w, T):
        """
        Also requires: environment temperature (density), so that buoyancy can be computed
        """
        g, gamma, mu = self.g, self.gamma, self.mu(r)

        B = (T-self.T_e(z))/T

        return 1./w*(g/(1+gamma)*B - mu*w**2)

    def dTdz__dQdz(self, z, r, w, T, q_v):
        g, cp_d, mu = self.g, self.cp_d, self.mu(r)

        Te_ = self.T_e(z)
        p_e = self.p_e(z)

        dTdz_s = -g/cp_d - mu*(T - Te_)

        # Intermediate values
        # TODO: entrainment of air
        T_s = T + dTdz_s
        qv__s = q_v

        dQ = self.microphysics(T=T_s, p=p_e, qv=qv__s)
        dq_v = dQ[0]

        dTdz__latent_heat = -self.Lv/cp_d*dq_v

        dTdz_ = dTdz_s + dTdz__latent_heat

        return dTdz_, dQ

    def drdz(self, z, r, w, T, dwdz_, dTdz_):
        g, R_d, mu = self.g, self.R_d, self.mu(r)

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
        self.beta = kwargs.pop('beta', 0.2)
        self.gamma = kwargs.pop('gamma', 0.5)
        self.D = kwargs.pop('D', 0.0)
        super(DryAirOnly, self).__init__(**kwargs)

    def mu(self, r):
        return self.beta/r

    def dwdz(self, z, r, w, T):
        """
        Also requires: environment temperature (density), so that buoyancy can be computed
        """
        g, gamma, mu, D = self.g, self.gamma, self.mu(r), self.D

        B = (T-self.T_e(z))/T

        return 1./w*(g/(1+gamma)*B - mu*w**2 - D*w**2/r)

    def dTdz(self, z, r, w, T):
        g, cp_d, mu = self.g, self.cp_d, self.mu(r)

        Te_ = self.T_e(z)

        return -g/cp_d - mu*(T - Te_)

    def drdz(self, z, r, w, T, dwdz_, dTdz_):
        g, R_d, mu = self.g, self.R_d, self.mu(r)

        return r/2. *( (g/(R_d*T) + 1./T*dTdz_) - 1./w * dwdz_ + mu/r)

    def dFdz(self, F, z):
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]

        dwdz_ = self.dwdz(z, r, w, T)
        dTdz_ = self.dTdz(z, r, w, T)
        drdz_ = self.drdz(z, r, w, T, dwdz_, dTdz_)

        return [drdz_, dwdz_, dTdz_, 0., 0., 0., 0.,]

class CCFM_v0(CloudModel):
    #import ccfm.version0

    def dFdz(self, F, z):
        raise NotImplemented

        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]
        q_v = F[Var.q_v]
        q_r = F[Var.q_r]
        q_l = F[Var.q_l]
        q_i = F[Var.q_i]

        p_e = self.p_e(z)

        warnings.warn("Environment assumed dry")
        T_e = self.T_e(z)
        qv_e = 0.0
        Tv_e = T_e*(1. + 0.61*qv_e)

        Tv_e = self.T_e(z)*(1. + 0.61*q_v)

        dwdz_ = ccfm.version0.ddv_z(vz_=w, q_l=q_l, q_r=q_r, q_i=q_i, q_s=0.0, t_ve=Tv_e, t_vc=Tv_c)
        dTdz_ = ccfm.version0.ddt_z(q_satw=q_satw, t_c=T, mu=mu, t_e=T_e, qv_e=qv_e, p_e=p_e, dpdz=dpdz)

        dpdz = 0.0
        warnings.warn("Should dpdz be zero?")
        drdz_ = ccfm.version0.ddr_z(ddt_z=dTdz_,ddv_z=dwdz_,t_ve=Tv_e,t_vc=Tv_c,rad=r,v_z=w)

        return [drdz_, dwdz_, dTdz_, 0., 0., 0., 0.,]
