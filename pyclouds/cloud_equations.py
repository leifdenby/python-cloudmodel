"""
Collection of cloud-model equations.
"""
import odespy
import warnings
import numpy as np
import plotting

from common import AttrDict, Var, default_constants
import cloud_microphysics


class CloudProfile():
    def __init__(self, F, z, cloud_model, extra_vars={}):
        self.F = F
        self.z = z
        self.cloud_model = cloud_model
        self.extra_vars = extra_vars

    def plot(self):
        return plotting.plot_profiles([self,], variables=['r', 'w', 'T',])

    def __str__(self):
        return str(self.cloud_model)


class CloudModel(object):
    ALLOWS_FREEZING = False

    def __init__(self, environment, constants=default_constants, microphysics=None):
        self.microphysics = microphysics
        self.environment = environment
        self.constants = AttrDict(constants)

        self.g = constants.get('g')
        self.R_d = constants.get('R_d')
        self.cp_d = constants.get('cp_d')

    def dFdz(F, z):
        raise NotImplemented

    def _stopping_criterion(self, F_solution, z, k):
        F_top = F_solution[k]
        z_top = z[k]

        if F_top[Var.T] > 300.0:
            print "Integration stopped: temperature got unphysically high"
            return True
        elif F_top[Var.w] < 0.0:
            print "Integration stopped: vertical velocity dropped to zero"
            print F_top[Var.w]
            return True
        elif F_top[Var.r] > 100e3:
            print "Integration stopped: cloud radius became unreasonably high (r>100km)"
            return True
        elif z_top > 30e3:
            print "Integration stopped: height reached too high (z > 30km)"
            return True
        elif np.any(np.isnan(F_top)):
            print "Integration stopped: solution became nan"
            if self.fail_on_nan:
                raise Exception("Solution became nan")
            return True
        else:
            return False

    def _stopping_criterion_time(self, F_solution, t, k):
        if hasattr(self, 'stop_integration'):
            return True

        F_top = F_solution[k]
        t_top = t[k]

        if F_top[Var.T] > 300.0:
            print "Integration stopped: temperature got unphysically high"
            return True
        elif F_top[Var.T] < 0.0:
            print "Integration stopped: temperature dropped below zero"
            return True
        elif F_top[Var.r] > 100e3:
            print "Integration stopped: cloud radius became unreasonably high (r>100km)"
            return True
        elif F_top[Var.z] < 0.0:
            print "Integration stopped: height below ground"
            return True
        elif F_top[Var.r] < 0.0:
            print "Integration stopped: radius dropped below zero"
            return True
        else:
            return False

    def dFdt(self, F, t):
        # print "--> ",
        # Var.print_formatted(F)
        w = F[Var.w]
        z = F[Var.z]

        try:
            dFdz = self.dFdz(F=F, z=z)
            dzdt = w
            dFdt = dFdz*dzdt
            dFdt[Var.z] = dzdt

            # print "dFdz",
            # Var.print_formatted(dFdz)

            # print "dFdt=",
            # Var.print_formatted(dFdt)
            # print
        except Exception as e:
            # TODO: stop integration when exception is raised
            self.stop_integration = True
            print "Error: %s" % e.message
            print "Exception at t=%fs" % t
            dFdt = np.zeros((Var.NUM,))

        return dFdt

    @staticmethod
    def _validate_initial_state(F):
        if F[Var.T] == 0.0:
            raise Exception("Should have temperature > 0.0")
        if F[Var.p] == 0.0:
            raise Exception("Should have pressure > 0.0")

    def integrate_in_time(self, initial_condition, t, SolverClass=odespy.RKFehlberg, stopping_criterion=None, tolerance=1e-3):
        warnings.warn("Integration in time is still experimental")

        self._validate_initial_state(initial_condition)

        solver = SolverClass(self.dFdt, rtol=0.0, atol=tolerance,)
        solver.set_initial_condition(np.array(initial_condition))

        if stopping_criterion is None:
            stopping_criterion=self._stopping_criterion_time

        F, t = solver.solve(t, stopping_criterion)

        return CloudProfile(F=F, z=F[:,Var.z], cloud_model=self)

    def integrate(self, initial_condition, z, SolverClass=odespy.RKFehlberg, stopping_criterion=None, tolerance=1e-3, fail_on_nan=False):
        self._validate_initial_state(initial_condition)

        rtol = 0.01
        atol_q = 1.0e-10
        min_step = 0.05
        atol = Var.make_state(p=100., T=0.1, q_v=atol_q, q_l=atol_q, w=0.1, r=10., q_r=atol_q, z=10., q_i=atol_q)

        solver = SolverClass(self.dFdz, rtol=rtol, atol=atol, min_step=min_step)
        solver.set_initial_condition(initial_condition)
        self.fail_on_nan = fail_on_nan

        mphys_extra_vars = {}
        self.microphysics.extra_vars = mphys_extra_vars

        if stopping_criterion is None:
            stopping_criterion=self._stopping_criterion
        F, z = solver.solve(z, stopping_criterion,)

        if self._stopping_criterion(F_solution=F, z=z, k=len(z)-1):
            F = F[:-1]
            z = z[:-1]

            for k, v in mphys_extra_vars.items():
                try:
                    mphys_extra_vars[k] = v[:-1]
                except IndexError:
                    pass

        return CloudProfile(F=F, z=z, cloud_model=self, extra_vars=mphys_extra_vars)


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

        T_e = self.environment.temp(z)

        B = (T-T_e)/T

        return 1./w*(g/(1+gamma)*B - mu*w**2 - D*w**2/r)

    def dTdz(self, z, r, w, T):
        g, cp_d, mu = self.g, self.cp_d, self.mu(r)

        Te_ = self.environment.temp(z)


        return -g/cp_d - mu*(T - Te_)

    def drdz(self, z, r, w, T, dwdz_, dTdz_):
        g, R_d, mu = self.g, self.R_d, self.mu(r)


        return r/2. *( (g/(R_d*T) + 1./T*dTdz_) - 1./w * dwdz_ + mu)

    def dFdz(self, F, z):
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]

        dwdz_ = self.dwdz(z, r, w, T)
        dTdz_ = self.dTdz(z, r, w, T)
        drdz_ = self.drdz(z, r, w, T, dwdz_, dTdz_)

        dFdz_ = np.zeros((Var.NUM,))
        dFdz_[Var.w] = dwdz_
        dFdz_[Var.T] = dTdz_
        dFdz_[Var.r] = drdz_

        return dFdz_

    def __str__(self):
        return r"DryAirEqns ($D=%g$, $\beta=%g$)" % (self.D, self.beta,)

class CCFM_v0(CloudModel):
    def dFdz(self, F, z):
        #import ccfm.version0

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
    def __init__(self, environment, gamma=0.5, D=0.506, beta=0.2, **kwargs):
        """
        gamma: virtual mass coefficient
        D: drag coefficient
        beta: entrainment coefficient

        Default values from Simpson & Wiggert 1969
        """
        self.entrain_moist_static_energy = kwargs.pop('entrain_moist_static_energy', True)
        super(FullThermodynamicsCloudEquations, self).__init__(environment=environment, **kwargs)

        self.gamma = gamma
        self.D = D
        self.beta = beta


    def cloud_mixture_density(self, p, T_c, qd_c, qv_c, ql_c, qr_c, qi_c):
        """
        Compute the mixture density from the import full equation of state

        Constants:
            R_d: specific gas constant of dry air
            R_v: specific gas constant of water vapour
            rho_l: density of liquid water
            rho_i: density of ice
        """
        R_d = self.constants.R_d
        R_v = self.constants.R_v
        rho_l = self.constants.rho_l
        rho_i = self.constants.rho_i

        rho_inv = (qd_c*R_d + qv_c*R_v)*T_c/p + (ql_c + qr_c)/rho_l + qi_c/rho_i
        
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
        R_d = self.constants.R_d
        R_v = self.constants.R_v

        R_s = (R_v*qv_c + R_d*qd_c)/(qd_c + qv_c)

        return p/(T_c*R_s)


    def dw_dz(self, p, w_c, r_c, T_c, qd_c, qv_c, ql_c, qi_c, qr_c, rho_e):
        """
        Momentum equation

        State variables:
            w: verticaly velocity
            r: cloud radius
            rho_c: cloud density
        """
        rho_c = self.cloud_mixture_density(p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c, ql_c=ql_c, qi_c=qi_c, qr_c=qr_c)

        g = self.constants.g

        B = (rho_e - rho_c)/rho_e

        mu = self.beta/r_c

        return 1.0/w_c * (g/(1.+self.gamma)*B - mu*w_c**2. - self.D*w_c**2./r_c)

    def dT_dz(self, r_c, T_c, qd_c, qv_c, ql_c, qi_c, dql_c__dz, dqi_c__dz, T_e, qv_e):
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
        cp_d = self.constants.cp_d
        cp_v = self.constants.cp_v
        L_v = self.constants.L_v
        L_s = self.constants.L_s
        g = self.constants.g

        if self.entrain_moist_static_energy:
            # heat capacity of cloud mixture
            c_cm_p = cp_d*qd_c + cp_v*(qv_c + ql_c + qi_c)

            qd_e = 1.0 - qv_e
            ql_e, qi_e = 0.0, 0.0

            c_em_p = cp_d*qd_e + cp_v*(qv_e + ql_e + qi_e)

            # difference in environment and cloud moist static energy
            Ds = (c_em_p*T_e + ql_e*L_v)\
                -(c_cm_p*T_c + ql_c*L_v)

            # XXX: There appears to be something wrong with the formulation that
            # includes liquid and ice, use just the liquid formulation for now
            # Ds = (c_em_p*T_e - ql_e*L_v - qi_e*L_s)\
                # -(c_cm_p*T_c - ql_c*L_v - qi_c*L_s)
        else:
            # previously the moist static energy was used for the cloud mixture
            # and for the environment only dry air was included which lead to
            # inconsistencies
            Ds = cp_d*T_e - cp_d*T_c

            c_cm_p = cp_d
            c_em_p = cp_d


        mu = self.beta/r_c

        return  -g/c_cm_p + mu*Ds/c_cm_p + L_v/c_cm_p*dql_c__dz + L_s/c_cm_p*dqi_c__dz

    def dr_dz(self, p, w_c, r_c, T_c, qd_c, qv_c, ql_c, qr_c, qi_c, dql_c__dz, dqi_c__dz, dTc_dz, dw_dz):
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
        R_d = self.constants.R_d
        R_v = self.constants.R_v
        g = self.constants.g
        rho_i = self.constants.rho_i
        rho_l = self.constants.rho_l

        # XXX: for now assume that specific concentration of in-cloud dry does not change
        dqd_c__dz = 0.0
        
        # XXX: assume that changes in water vapour are exactly opposite to the increase in liquid water and ice
        dqv_c__dz = - dql_c__dz - dqi_c__dz

        # in-cloud mixture density
        rho_c = self.cloud_mixture_density(p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c, ql_c=ql_c, qi_c=qi_c, qr_c=qr_c)

        # in-cloud gas density
        rho_cg = self._cloud_gas_density_from_eos(p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c)

        # effective gas constant
        Rs_c = (R_v*qv_c + R_d*qd_c)/(qv_c + qd_c)

        # in-cloud specific constant of gas constituents
        qg_c = qd_c + qv_c

        mu = self.beta/r_c


        return r_c/2.0*(\
                        qg_c*rho_c/rho_cg * rho_c/rho_cg * g/(Rs_c*T_c)\
                        + qg_c*rho_c/rho_cg * 1./T_c*dTc_dz\
                        + rho_c/(rho_cg*Rs_c) * (dqv_c__dz*R_v + dqd_c__dz*R_d)\
                        + rho_c/rho_i*dqi_c__dz\
                        + rho_c/rho_l*dql_c__dz\
                        + mu #  1./M*dM_dz\
                        - 1./w_c*dw_dz\
                      )

    def dFdz(self, F, z):
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]
        q_v = F[Var.q_v]
        q_l = F[Var.q_l]
        q_i = F[Var.q_i]
        q_r = F[Var.q_r]
        q_d = 1.0 - q_v - q_l - q_i - q_r

        # cloud is assumed to be at same pressure as in environment
        p = self.environment.p(z)
        # NB: need to make sure that pressure is set since the microphysics
        # needs it and we don't integrate a pressure gradient, instead we
        # assume pressure balance with environment
        F[Var.p] = p

        # print "z=", z,
        # Var.print_formatted(F)

        rho_e = self.environment.rho(z)
        T_e = self.environment.temp(z)

        dFdz_ = np.zeros((Var.NUM,))

        # 1. Estimate change in vertical velocity with initial state
        dwdz_ = self.dw_dz(p=p, w_c=w, r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qr_c=q_r, qi_c=q_i, rho_e=rho_e)

        dFdz_[Var.w] = dwdz_

        # 2. Estimate temperature change forgetting about phase-changes for now (i.e. considering only adiabatic adjustment and entrainment)
        try:
            qv_e = self.environment.rel_humidity(z)*self.microphysics.parameterisations.pv_sat.qv_sat(T=T_e, p=p)
        except AttributeError:
            warnings.warn("It seems the environmental profile doesn't define a relative humidity so we'll assume it's dry")
            qv_e = 0.0

        dTdz_s = self.dT_dz(r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i, dql_c__dz=0.0, dqi_c__dz=0.0,
                            T_e=T_e, qv_e=qv_e)

        dFdz_[Var.T] = dTdz_s

        # 3. estimate new state from phase changes predicted by microphysics
        F[Var.p] = p  # make sure that pressure is set since the microphysics needs it
        dFdt_micro = self.microphysics.dFdt(F, t=z) # passing in `z` as `t` doesn't make sense physically, but it will just be used for plotting later

        dFdz_micro = dFdt_micro/w  # w = dz/dt
        dFdz_ += dFdz_micro

        dql_c__dz = dFdz_micro[Var.q_l]
        dqi_c__dz = dFdz_micro[Var.q_i]
        dTdz_ = dTdz_s + dFdz_micro[Var.T]

        dFdz_[Var.T] = dTdz_

        # 4. Use post microphysics state (and phase changes from microphysics) to estimate radius change
        drdz_ = self.dr_dz(p=p, w_c=w, r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qr_c=q_r,
                           qi_c=q_i, dql_c__dz=dql_c__dz, dqi_c__dz=dqi_c__dz, dTc_dz=dTdz_, dw_dz=dwdz_)

        dFdz_[Var.r] = drdz_

        # calculate adiabatic lapse rate from equation on wikipedia
        # from pyclouds.common import default_constants

        # pv_sat = self.microphysics.parameterisations.pv_sat(T=T)
        # g = default_constants.get('g')
        # L_v = default_constants.get('L_v')
        # R_v = default_constants.get('R_v')
        # R_d = default_constants.get('R_d')
        # eps = R_d/R_v
        # cp_d = default_constants.get('cp_d')
        
        # r = eps*pv_sat/(p - pv_sat)

        # dTdz_moist = g*(1. + L_v*r/(R_d*T))/(cp_d + L_v**2.*r*eps/(R_d*T**2.))

        # qv_sat = self.microphysics.parameterisations.pv_sat.qv_sat(T=T, p=p)
        # print "2> qv={}, qv_sat={}, T={}, p={}".format(q_v, qv_sat, T, p)
        # print z, dTdz_s, dFdz_micro[Var.T], dTdz_moist, 'Sw=', q_v/qv_sat
        # print "ddz", 
        # Var.print_formatted(dFdz_)
        # print "dz_max=",
        # Var.print_formatted(F/dFdz_)
        # print ">> dF (dz=5m)=",
        # Var.print_formatted(dFdz_*5.)
        # Var.print_formatted(F + dFdz_*5.)
        # print

        # raw_input()

        return dFdz_

    def __str__(self):
        return r"FullSpecConcEqns ($D=%g$, $\beta=%g$), mu-phys: %s" % (self.D, self.beta, str(self.microphysics))

class FullEquationsSatMicrophysics(FullThermodynamicsCloudEquations):
    """
    Same as the full model, but using the moist-adjustment microphysics is hardcoded in
    """
    def __init__(self, *args, **kwargs):
        if 'microphysics' in kwargs:
            warnings.warn("%s is hardcoded to use the moist adjustment microphysics")
            kwargs['microphysics'] = cloud_microphysics.MoistAdjustmentMicrophysics()
        super(FullEquationsSatMicrophysics, self).__init__(*args, **kwargs)


    def dFdz(self, F, z):
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]
        q_v = F[Var.q_v]
        q_l = F[Var.q_l]
        q_i = F[Var.q_i]
        q_r = F[Var.q_r]
        q_d = 1.0 - q_v - q_l - q_i - q_r

        if q_r != 0.0:
            raise NotImplementedError

        # cloud is assumed to be at same pressure as in environment
        p = self.environment.p(z)

        dFdz_ = np.zeros((Var.NUM,))
        

        # 1. Estimate change in vertical velocity with initial state
        rho_e = self.environment.rho(z)
        dwdz_ = self.dw_dz(p=p, w_c=w, r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i, rho_e=rho_e)

        dFdz_[Var.w] = dwdz_

        # 2. Estimate temperature change forgetting about phase-changes for now (i.e. considering only adiabatic adjustment and entrainment)
        T_e = self.environment.temp(z)
        dTdz_s = self.dT_dz(r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i, dql_c__dz=0.0, dqi_c__dz=0.0, T_e=T_e)

        F_s = np.copy(F)
        F_s[Var.T] += dTdz_s


        # 3. With new temperature estimate new state from phase changes predicted by microphysics
        F_s = self.microphysics._calc_adjusted_state(F_s, iterations=2)

        dTdz_ = F_s[Var.T] - F[Var.T]
        dFdz_[Var.T] = dTdz_

        dql_c__dz = F_s[Var.q_l] - F[Var.q_l]
        dqi_c__dz = F_s[Var.q_i] - F[Var.q_i]
        dFdz_[Var.q_v] = - dql_c__dz - dqi_c__dz
        dFdz_[Var.q_l] = dql_c__dz
        dFdz_[Var.q_i] = dqi_c__dz

        # 4. Use post microphysics state (and phase changes from microphysics) to estimate radius change
        drdz_ = self.dr_dz(p=p, w_c=w, r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i, dql_c__dz=dql_c__dz, dqi_c__dz=dqi_c__dz, dTc_dz=dTdz_, dw_dz=dwdz_)

        dFdz_[Var.r] = drdz_

        return dFdz_

class FixedRiseRateParcel(CloudModel):
    """
    Simple model to test the lapse rate of a parcel that as is forced to rise
    with a constant velocity
    """
    def __init__(self, w0, beta, **kwargs):
        """
        gamma: virtual mass coefficient
        D: drag coefficient
        beta: entrainment coefficient
        """
        super(self.__class__, self).__init__(**kwargs)

        self.w0 = w0
        self.beta = beta

    def _stopping_criterion(self, F_solution, z, k):
        F_top = F_solution[k]
        z_top = z[k]

        if z_top > 5000.0:
            print "Integration stopped: max altitude reached"
            return True
        else:
            return False

    def cloud_mixture_density(self, p, T_c, qd_c, qv_c, ql_c, qr_c, qi_c):
        """
        Compute the mixture density from the import full equation of state

        Constants:
            R_d: specific gas constant of dry air
            R_v: specific gas constant of water vapour
            rho_l: density of liquid water
            rho_i: density of ice
        """
        R_d = self.constants.R_d
        R_v = self.constants.R_v
        rho_l = self.constants.rho_l
        rho_i = self.constants.rho_i

        rho_inv = (qd_c*R_d + qv_c*R_v)*T_c/p + (ql_c + qr_c)/rho_l + qi_c/rho_i
        
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
        R_d = self.constants.R_d
        R_v = self.constants.R_v

        R_s = (R_v*qv_c + R_d*qd_c)/(qd_c + qv_c)

        return p/(T_c*R_s)


    def dT_dz(self, r_c, T_c, qd_c, qv_c, ql_c, qi_c, dql_c__dz, dqi_c__dz, T_e):
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
        cp_d = self.constants.cp_d
        cp_v = self.constants.cp_v
        L_v = self.constants.L_v
        L_s = self.constants.L_s
        g = self.constants.g

        # heat capacity of cloud mixture
        c_cm_p = cp_d*qd_c + cp_v*(qv_c + ql_c + qi_c)

        # heat capacity of environment mixture
        # qd_e, qv_e, ql_e, qi_e = self.qd_e, self.qv_e, self.ql_e, self.qi_e
        # TODO: Allow definition of moist environment
        qd_e, qv_e, ql_e, qi_e = 1.0, 0.0, 0.0, 0.0

        c_em_p = cp_d*qd_e + cp_v*(qv_e + ql_e + qi_e)

        # difference in environment and cloud moist static energy
        # XXX: There appears to be something wrong with the formulation that
        # includes liquid and ice, use just the liquid formulation for now
        # Ds = (c_em_p*T_e - ql_e*L_v - qi_e*L_s)\
            # -(c_cm_p*T_c - ql_c*L_v - qi_c*L_s)
        Ds = (c_em_p*T_e + ql_e*L_v)\
            -(c_cm_p*T_c + ql_c*L_v)

        mu = self.beta/r_c


        return  -g/c_cm_p + mu*Ds/c_cm_p + L_v/c_cm_p*dql_c__dz + L_s/c_cm_p*dqi_c__dz


    def dFdz(self, F, z):
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]
        q_v = F[Var.q_v]
        q_l = F[Var.q_l]
        q_i = F[Var.q_i]
        q_r = F[Var.q_r]
        q_d = 1.0 - q_v - q_l - q_i - q_r

        # cloud is assumed to be at same pressure as in environment
        p = self.environment.p(z)

        rho_e = self.environment.rho(z)
        T_e = self.environment.temp(z)

        dFdz_ = np.zeros((Var.NUM,))

        dFdz_[Var.w] = 0.0

        # 2. Estimate temperature change forgetting about phase-changes for now (i.e. considering only adiabatic adjustment and entrainment)
        dTdz_s = self.dT_dz(r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i, dql_c__dz=0.0, dqi_c__dz=0.0, T_e=T_e)

        dFdz_[Var.T] = dTdz_s

        # 3. estimate new state from phase changes predicted by microphysics
        dFdt_micro = self.microphysics.dFdt(F, t=None)
        dFdz_micro = dFdt_micro/w  # w = dz/dt
        dFdz_ += dFdz_micro

        dFdz_[Var.r] = 0.0

        return dFdz_

    def __str__(self):
        return r"Fixed vertical velocity: w0={w0}m/s, mu-phys: {mphys}".format(w0=self.w0, mphys=self.microphysics)
