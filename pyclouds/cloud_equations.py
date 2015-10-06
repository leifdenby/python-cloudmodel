"""
Collection of cloud-model equations.
"""
import odespy
import warnings

# r, w, T, q_v, q_r, q_l, q_i
class Var:
    p = 0
    r = 1
    w = 2
    T = 3
    q_v = 4
    q_r = 5
    q_l = 6
    q_i = 7

class VarQ:
    q_v = 0
    q_r = 1
    q_l = 2
    q_i = 3


class CloudModel(object):
    ALLOWS_FREEZING = False

    def __init__(self, environment, constants, microphysics=None):
        self.microphysics = microphysics

        self.T_e = environment.get('T_e')
        self.p_e = environment.get('p_e')
        self.g = constants.get('g')
        self.R_d = constants.get('R_d')
        self.cp_d = constants.get('cp_d')

    def dFdz(F, z):
        raise NotImplemented

    def _stopping_criterion(self, F_solution, z, k):
        F_top = F_solution[k]

        if not self.ALLOWS_FREEZING and F_top[Var.T] < 0.0:
            print "Integration stopped: model does not allow freezing"
            return True
        elif F_top[Var.T] > 300.0:
            print "Integration stopped: temperature got unphysically high"
            return True
        elif F_top[Var.w] < 0.0:
            print "Integration stopped: vertical velocity dropped to zero"
            return True
        else:
            return False

    def integrate(self, initial_condition, z, SolverClass=odespy.RKFehlberg):
        solver = SolverClass(self.dFdz, rtol=0.0, atol=1e-6,)
        solver.set_initial_condition(initial_condition)

        F, z = solver.solve(z, self._stopping_criterion)

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
        dq_v = dQ[VarQ.q_v]

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

        # p, r, w, T, q_v, q_r, q_l, q_i
        return [0.0, drdz_, dwdz_, dTdz_, dQdz_[0], dQdz_[1], 0., 0.,]

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


class FullThermodynamicsCloudEquations(CloudModel):
    """
    This class represents the full cloud equations where constituents written
    in terms of specific concentrations

    The model assumes that:
    1. The cloud temperature at given height is the same as the environment at that height
    2. The pressure in-cloud and in-environment is in hydrostatic equilibrium
    ...

    """
    def __init__(self, gamma, D, beta, rho0, **kwargs):
        """
        gamma: virtual mass coefficient
        D: drag coefficient
        beta: entrainment coefficient
        """
        super(self, FullThermodynamicsCloudEquations).__init__(**kwargs)

        self.gamma = gamma
        self.D = D
        self.beta = beta

    def _cloud_mixture_density_from_eos(self, p, T_c, qd_c, qv_c, ql_c, qi_c):
        """
        Compute the mixture density from the import full equation of state

        Constants:
            R_d: specific gas constant of dry air
            R_v: specific gas constant of water vapour
            rho_l: density of liquid water
            rho_i: density of ice
        """
        rho_inv = (qd_c*R_d + qv_c*R_v)*T_c/p + ql_c/rho_l + qi_c/rho_i
        
        return 1.0/rho_inv

    def _cloud_gas_density_from_eos(self, p, T_c, qd_c, qv_c):
        """
        Compute the gas density from equation of state. This is the same as the
        full equation of state, but rearranged to have the form of a
        traditional ideal gas equation of state.

        Constants:
            R_d: specific gas constant of dry air
            R_v: specific gas constant of water vapour
            R_s: effective specific gas constant of gas mixture
        """
        R_s = (R_v*qv_c + R_d*qd_c)/(qd_c + qv_c)

        return p/(T*R_s)

    def dp_dz(self, p, T_c, qd_c, qv_c, ql_c, qi_c):
        rho_c = self._cloud_mixture_density_from_eos(p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c, ql_c=ql_c, qi_c=qi_c)

        return -rho_c*g

    def dw_dz(self, p, w_c, r_c, T_c, qd_c, qv_c, ql_c, qi_c):
        """
        Momentum equation

        State variables:
            w: verticaly velocity
            r: cloud radius
            rho_c: cloud density
        """
        rho_c = self._cloud_mixture_density_from_eos(p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c, ql_c=ql_c, qi_c=qi_c)

        return g/(1.+self.gamma)(rho_c - self.rho0)/rho_c - self.beta/r_c*w_c**2. - self.D*w_c**2./r_c

    def dT_dz(self, p, w_c, r_c, T_c, qd_c, qv_c, ql_c, qi_c, dql_c__dz, dqi_c__dz):
        """
        Constants:
            cp_d: heat capacity of dry air at constant pressure
            cp_v: heat capacity of liquid water at constant pressure
            L_v: latent heat of vapourisation (vapour -> liquid)
            L_s: latent heat sublimation (vapour -> solid)
            
        State variables:
            qd_c, qv_c, ql_c, qi_c: in-cloud dry air, water vapour, liquid water, ice
            T_c: in-cloud absolute temp

            dql_c__dz, qdi_c__dz: vertical in

            qd_c, qv_c, ql_c, qi_c: environment (constant) dry air, water vapour, liquid water, ice
            T_e: environment absolute temp
        """
        # heat capacity of cloud mixture
        c_cm_p = cp_d*qd_c + cv_d*(qv_c + ql_c + qi_c)

        # heat capacity of environment mixture
        qd_e, qv_e, ql_e, qi_e = self.qd_e, self.qv_e, self.ql_e, self.qi_e
        c_em_p = cp_d*qd_e + cv_d*(qv_e + ql_e + qi_e)

        # difference in environment and cloud moist static energy
        Ds = (c_em_p*T_e - ql_e*L_v - qi_e*L_s)\
            -(c_cm_p*T_c - ql_c*L_v - qi_c*L_s)


        return  g/c_cm_p + self.beta/r*Ds + Lv/c_cm_p*dql_c__dz + Ls/c_cm_p*dqi_c__dz

    def dr_dz(self, p, w_c, r_c, T_c, qd_c, qv_c, ql_c, qi_c, dql_c__dz, dqi_c__dz):
        """
        Mass conservation equation

        Constants
            R_d: specific heat capacity of dry air
            R_v: specific heat capacity of water vapour


        State variables:
            qc_g: specific concentration of gas species in-cloud
            T_c: absolute temperature in cloud
            rho_c: density in-cloud
            rho_cg: density of gas mixture (wrt gas volume)

            (assumed constant)
            rho_i: in-cloud density of ice
            rho_l: in-cloud density of liquid water
        """
        # XXX: for now assume that specific concentration of in-cloud dry does not change
        dqd_c__dz = 0.0
        
        # XXX: assume that changes in water vapour are exactly opposite to the increase in liquid water and ice
        dqv_c__dz = - dql_c__dz - dqi_c__dz

        # in-cloud mixture density
        rho_c = self._cloud_mixture_density_from_eos(p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c, ql_c=ql_c, qi_c=qi_c)

        # in-cloud gas density
        rho_cg = self._cloud_gas_density_from_eos(p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c)

        # effective gas constant
        Rs_c = (R_v*qv_c + R_d*qd_c)/(qv_c + qd_c)

        # in-cloud specific constant of gas constituents
        qg_c = qd_c + qv_c

        return r/2.0*(\
                        qg_c*rho_c/rho_g * rho_c/rho_cg * g/Rs_c*T_c\
                        + qg_c*rho_c/rho_g * 1./T_c*dTc_dz\
                        + rho_c/(rho_cg*Rs_c) * (dqv_c__dz*Rv + dqd_c__dz*Rd)\
                        + rho_c/rho_i*dqc_i__dz\
                        + rho_c/rho_l*dqc_l__dz\
                        + self.beta/r #  1./M*dM_dz\
                        - 1./w*dw_dz\
                      )

    def dFdz(self, F, z):
        p = F[Var.p]
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]
        q_v = F[Var.q_v]
        q_l = F[Var.q_l]
        q_i = F[Var.q_i]
        q_r = F[Var.q_r]
        q_d = 1.0 - q_v - q_l - q_i - q_r

        dFdz_ = np.zeros(F)
        

        # 1. Estimate change in vertical velocity with initial state
        dwdz_ = self.dw_dz(p=p, w_c=w, r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i)

        dFdz_[Var.w] = dFdz_

        # 2. Estimate temperature change forgetting about phase-changes for now (i.e. considering only adiabatic adjustment and entrainment)
        dTdz_s = self.dT_dz(p=p, w_c=w, r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i, dql_c__dz=0.0, dqi_c__dz=0.0)

        F_s = np.copy(F)
        F_s[Var.T] += dTdz_s


        # 3. With new temperature estimate new state from phase changes predicted by microphysics
        F_s = self.microphysics(F_s)

        dTdz_ = F[Var.T] - F_s[Var.T]
        dFdz_[Var.T] = dTdz_

        dql_c__dz = F[Var.q_l] - F_t[Var.q_l]
        dqi_c__dz = F[Var.q_i] - F_t[Var.q_i]
        dFdz_[Var.d_v] = - dql_c__dz - dqi_c__dz
        dFdz_[Var.d_l] = dql_c__dz
        dFdz_[Var.d_i] = dqi_c__dz

        # 4. Use post microphysics state (and phase changes from microphysics) to estimate radius change
        drdz_ = self.dr_dz(p=p, w_c=w, r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i, dql_c__dz=dql_c__dz, dqi_c__dz=dqi_c__dz)

        dFdz_[Var.r] = drdz_

        # 5. Estimate pressure change from original state
        # XXX: Shouldn't we use this new temperature in the steps above? Or at least use the updated state variables to compute the new pressure?
        dpdz_ = self.dp_dz(p=p, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i)

        dFdz_[Var.p] = dpdz_

        return dFdz_
