"""
Collection of cloud-model equations.
"""
import warnings
import numpy as np
from scipy.constants import pi

from .. import AttrDict, Var
from ..reference.constants import default_constants
from . import microphysics as cloud_microphysics
from ..integration.parcel import ParcelModelIntegrator


class CloudModel(object):
    ALLOWS_FREEZING = False

    def __init__(self, environment, constants=default_constants, microphysics=None):
        self.microphysics = microphysics
        self.environment = environment
        self.constants = AttrDict(constants)

        self.g = constants.get("g")
        self.R_d = constants.get("R_d")
        self.cp_d = constants.get("cp_d")

        self.extra_vars = {}

    def integrate(self, *args, **kwargs):
        integrator = ParcelModelIntegrator(cloud_model=self)
        return integrator(*args, **kwargs)

    def dFdz(z, F):
        raise NotImplemented

    def dFdt(self, t, F):
        w = F[Var.w]
        z = F[Var.z]

        try:
            dFdz = self.dFdz(F=F, z=z)
            dzdt = w
            dFdt = dFdz * dzdt
            dFdt[Var.z] = dzdt
        except Exception as e:
            # TODO: stop integration when exception is raised
            self.stop_integration = True
            print("Error: %s" % e.message)
            print("Exception at t=%fs" % t)
            dFdt = np.zeros((Var.NUM,))

        return dFdt


class NoMicrophysicsNoEntrainment(CloudModel):
    mu = 0.0  # entrainment rate [1/m]
    gamma = 0.5  # virtual mass coefficient [1]
    D = 0.0  # drag coefficient [1]
    Lv = 2.5008e6  # latent heat of vapourisation [J/kg]

    def dwdz(self, z, r, w, T, Te):
        """
        Also requires: environment temperature (density), so that buoyancy can be computed
        """
        g = self.constants.g
        gamma = self.gamma
        mu = self.mu
        D = self.D

        B = (T - Te) / T

        return 1.0 / w * (g / (1 + gamma) * B - mu * w ** 2 - D * w ** 2 / r)

    def dTdz(self, z, r, w, T, Te, p):
        g = self.constants.g
        mu = self.mu
        Lv = self.Lv
        R_d = self.constants.R_d
        cp_d = self.constants.cp_d

        qsat_w = self.microphysics.parameterisations.pv_sat.qv_sat(T=T, p=p)
        return -g / cp_d * (1 + (Lv * qsat_w) / (R_d * T)) - mu * (Te - T)

    def drdz(self, z, r, w, T, dwdz_, dTdz_):
        # dwdz_ = dwdz(z, r, w, T)
        # dTdz_ = dTdz(z, r, w, T)
        g = self.constants.g
        R_d = self.constants.R_d
        mu = self.mu

        return 2.0 / r * ((g / (R_d * T) + 1.0 / T * dTdz_) - 1.0 / w * dwdz_ + mu / r)

    def dFdz(self, F, z):
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]

        Te = self.environment.temp(z)

        dwdz_ = self.dwdz(z, r, w, T, Te)
        dTdz_ = self.dTdz(z, r, w, T, Te)
        drdz_ = self.drdz(z, r, w, T, dwdz_, dTdz_)

        return [
            drdz_,
            dwdz_,
            dTdz_,
            0.0,
            0.0,
            0.0,
            0.0,
        ]


class Wagner2009(CloudModel):
    """
    Moist saturation cloud-model, from Till's PhD thesis
    """

    Lv = 2.5008e6  # latent heat of vapourisation [J/kg]
    beta = 0.2  # entraiment-rate coefficient [1]
    gamma = 0.5  # virtual-mass coefficient [1]

    def __init__(self, **kwargs):
        self.microphysics = kwargs.pop("microphysics")
        self.qv_e = kwargs["environment"].get("qv_e")
        super(Wagner2009, self).__init__(**kwargs)

    def mu(self, r):
        return self.beta / r

    def dwdz(self, z, r, w, T):
        """
        Also requires: environment temperature (density), so that buoyancy can be computed
        """
        g, gamma, mu = self.g, self.gamma, self.mu(r)

        B = (T - self.T_e(z)) / T

        return 1.0 / w * (g / (1 + gamma) * B - mu * w ** 2)

    def dTdz__dQdz(self, z, r, w, T, q_v):
        g, cp_d, mu = self.g, self.cp_d, self.mu(r)

        Te_ = self.T_e(z)
        p_e = self.p_e(z)

        dTdz_s = -g / cp_d - mu * (T - Te_)

        # Intermediate values
        # TODO: entrainment of air
        T_s = T + dTdz_s
        qv__s = q_v

        dQ = self.microphysics(T=T_s, p=p_e, qv=qv__s)
        dq_v = dQ[Var.q_v]

        dTdz__latent_heat = -self.Lv / cp_d * dq_v

        dTdz_ = dTdz_s + dTdz__latent_heat

        return dTdz_, dQ

    def drdz(self, z, r, w, T, dwdz_, dTdz_):
        g, R_d, mu = self.g, self.R_d, self.mu(r)

        return r / 2.0 * ((g / (R_d * T) + 1.0 / T * dTdz_) - 1.0 / w * dwdz_ + mu / r)

    def dFdz(self, z, F):
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]
        q_v = F[Var.q_v]

        dwdz_ = self.self.dwdz(z, r, w, T)
        dTdz_, dQdz_ = self.dTdz__dQdz(z, r, w, T, q_v)
        drdz_ = self.drdz(z, r, w, T, dwdz_, dTdz_)

        # p, r, w, T, q_v, q_r, q_l, q_i
        return [
            0.0,
            drdz_,
            dwdz_,
            dTdz_,
            dQdz_[0],
            dQdz_[1],
            0.0,
            0.0,
        ]


class DryAirOnly(CloudModel):
    def __init__(self, **kwargs):
        """
        mu: entrainment rate [1/m]
        gamma: virtual mass coefficient [1]
        C_D: drag coefficient [1]
        """
        self.beta = kwargs.pop("beta", 0.2)
        self.gamma = kwargs.pop("gamma", 0.5)
        self.C_D = kwargs.pop("C_D", 0.0)

        super(DryAirOnly, self).__init__(**kwargs)

    def mu(self, r):
        return self.beta / r

    def dwdz(self, z, r, w, T):
        """
        Also requires: environment temperature (density), so that buoyancy can be computed
        """
        g, gamma, mu, C_D = self.g, self.gamma, self.mu(r), self.C_D

        T_e = self.environment.temp(z)

        B = (T - T_e) / T

        return 1.0 / w * (g / (1 + gamma) * B - mu * w ** 2 - C_D * w ** 2 / r)

    def dTdz(self, z, r, w, T):
        g, cp_d, mu = self.g, self.cp_d, self.mu(r)

        Te_ = self.environment.temp(z)

        return -g / cp_d - mu * (T - Te_)

    def drdz(self, z, r, w, T, dwdz_, dTdz_):
        g, R_d, mu = self.g, self.R_d, self.mu(r)

        return r / 2.0 * ((g / (R_d * T) + 1.0 / T * dTdz_) - 1.0 / w * dwdz_ + mu)

    def dFdz(self, z, F):
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
        return r"DryAirEqns ($C_D=%g$, $\beta=%g$)" % (
            self.C_D,
            self.beta,
        )


class FullThermodynamicsCloudEquations(CloudModel):
    """
    This class represents the full cloud equations where constituents written
    in terms of specific concentrations

    The model assumes that:
    1. The cloud temperature at given height is the same as the environment at that height
    2. The pressure in-cloud and in-environment is in hydrostatic equilibrium
    ...

    """

    def __init__(
        self,
        environment,
        microphysics,
        gamma=0.5,
        C_D=0.506,
        beta=0.2,
        l_pr=100.0,
        **kwargs
    ):
        """
        gamma: virtual mass coefficient
        D: drag coefficient
        beta: entrainment coefficient

        Default values from Simpson & Wiggert 1969

        l_pr: rain-out vertical length-scale
        """
        self.entrain_moist_static_energy = kwargs.pop(
            "entrain_moist_static_energy", True
        )
        self.entrain_liquid_static_energy = kwargs.pop(
            "entrain_liquid_static_energy", False
        )
        self.entrain_hydrometeors = kwargs.pop("entrain_hydrometeors", True)
        self.temperature_dependent_latent_heats = kwargs.pop(
            "temperature_dependent_latent_heats", True
        )
        super(FullThermodynamicsCloudEquations, self).__init__(
            environment=environment, microphysics=microphysics, **kwargs
        )

        if self.entrain_moist_static_energy and self.entrain_liquid_static_energy:
            raise Exception("Can only entrain either moist or liquid static energy")

        self.gamma = gamma
        self.C_D = C_D
        self.beta = beta
        self.l_pr = l_pr

    def _mu(self, r, w, B, z):
        return self.beta / r

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

        rho_inv = (
            (qd_c * R_d + qv_c * R_v) * T_c / p + (ql_c + qr_c) / rho_l + qi_c / rho_i
        )

        return 1.0 / rho_inv

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

        R_s = (R_v * qv_c + R_d * qd_c) / (qd_c + qv_c)

        return p / (T_c * R_s)

    def dw_dz(self, p, w_c, r_c, T_c, qd_c, qv_c, ql_c, qi_c, qr_c, rho_e, mu):
        """
        Momentum equation

        State variables:
            w: verticaly velocity
            r: cloud radius
            rho_c: cloud density
        """
        rho_c = self.cloud_mixture_density(
            p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c, ql_c=ql_c, qi_c=qi_c, qr_c=qr_c
        )

        g = self.constants.g

        B = (rho_e - rho_c) / rho_e

        return (
            1.0
            / w_c
            * (
                g / (1.0 + self.gamma) * B
                - mu * w_c ** 2.0
                - 3.0 / 8.0 * self.C_D * w_c ** 2.0 / r_c
            )
        )

    def dT_dz(
        self,
        r_c,
        T_c,
        qd_c,
        qv_c,
        ql_c,
        qi_c,
        dql_c__dz,
        dqi_c__dz,
        dqv_c__dz,
        T_e,
        qv_e,
        mu,
    ):
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
        cp_l = self.constants.cp_l
        cp_i = self.constants.cp_i

        g = self.constants.g

        L_v0 = self.constants.L_v
        L_f0 = self.constants.L_f
        L_s0 = self.constants.L_s

        if self.temperature_dependent_latent_heats:
            T00 = 273.15
            L_v = L_v0 + (cp_v - cp_l) * (T_c - T00)
            L_s = L_s0 + (cp_v - cp_i) * (T_c - T00)
            L_f = L_v0 + (cp_l - cp_i) * (T_c - T00)
        else:
            L_v = L_v0
            L_s = L_s0
            L_f = L_v0

        if self.entrain_moist_static_energy:
            if qi_c != 0.0 and not np.isnan(qi_c):
                raise NotImplementedError

            c_cm_p = cp_d * qd_c + cp_l * (qv_c + ql_c + qi_c)

            qd_e = 1.0 - qv_e
            ql_e, qi_e, qr_e = 0.0, 0.0, 0.0

            c_em_p = cp_d * qd_e + cp_l * (qv_e + ql_e + qi_e)

            # difference in environment and cloud moist static energy
            Ds = (c_em_p * T_e + qv_e * L_v - qi_e * L_f) - (
                c_cm_p * T_c + qv_c * L_v - qi_c * L_f
            )

            dqd_c__dz = -(dqv_c__dz + dql_c__dz + dqi_c__dz)

            # heat changes due to phase changes (and entrainment of species)
            dedz_q = (
                T_c * cp_d * dqd_c__dz
                + (T_c * cp_l + L_v) * dqv_c__dz
                + T_c * cp_l * dql_c__dz
                + (T_c * cp_i - L_f) * dqi_c__dz
            )

            # *actual* mixture heat capacity (if temperature dependency of
            # latent heats are considered in derivation then the *actual*
            # mixture heat capacity comes out)
            c_cm = cp_d * qd_c + qv_c * cp_v + ql_c * cp_l + qi_c * cp_i

            return -g / c_cm + mu * Ds / c_cm - 1.0 / c_cm * dedz_q
        elif self.entrain_liquid_static_energy:
            # XXX: This is *actually* liquid static energy
            # heat capacity of mixture with all moisture in the vapour phase
            qr_c = 0.0  # TODO: is this reasonable?
            c_cm_p = cp_d * qd_c + cp_v * (qv_c + ql_c + qr_c + qi_c)

            qd_e = 1.0 - qv_e
            ql_e, qr_e, qi_e = 0.0, 0.0, 0.0

            c_em_p = cp_d * qd_e + cp_v * (qv_e + ql_e + qr_e + qi_e)

            # difference in environment and cloud moist static energy
            Ds = (c_em_p * T_e - ql_e * L_v) - (c_cm_p * T_c - ql_c * L_v)

            # XXX: There appears to be something wrong with the formulation that
            # includes liquid and ice, use just the liquid formulation for now
            # Ds = (c_em_p*T_e - ql_e*L_v - qi_e*L_s)\
            # -(c_cm_p*T_c - ql_c*L_v - qi_c*L_s)

            return (
                -g / c_cm_p
                + mu * Ds / c_cm_p
                + L_v / c_cm_p * dql_c__dz
                + L_s / c_cm_p * dqi_c__dz
            )
        else:
            return -g / cp_d

    def dr_dz(
        self,
        p,
        w_c,
        r_c,
        T_c,
        qd_c,
        qv_c,
        ql_c,
        qr_c,
        qi_c,
        dqv_c__dz,
        dql_c__dz,
        dqi_c__dz,
        dTc_dz,
        dw_dz,
        mu,
    ):
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

        # Total specific concentration stays unchanged
        dqd_c__dz = -(dqv_c__dz + dql_c__dz + dqi_c__dz)

        # in-cloud mixture density
        rho_c = self.cloud_mixture_density(
            p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c, ql_c=ql_c, qi_c=qi_c, qr_c=qr_c
        )

        # in-cloud gas density
        rho_cg = self._cloud_gas_density_from_eos(p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c)

        # effective gas constant
        Rs_c = (R_v * qv_c + R_d * qd_c) / (qv_c + qd_c)

        # in-cloud specific constant of gas constituents
        qg_c = qd_c + qv_c

        return (
            r_c
            / 2.0
            * (
                qg_c * rho_c / rho_cg * rho_c / rho_cg * g / (Rs_c * T_c)
                + qg_c * rho_c / rho_cg * 1.0 / T_c * dTc_dz
                + rho_c / (rho_cg * Rs_c) * (dqv_c__dz * R_v + dqd_c__dz * R_d)
                + rho_c / rho_i * dqi_c__dz
                + rho_c / rho_l * dql_c__dz
                + mu  #  1./M*dM_dz\
                - 1.0 / w_c * dw_dz
            )
        )

    def dqr_dz__rainout(self, rho_c, q_r, w):
        """
        Estimate rate at which rain-droplets leave the cloudy air parcel. Based
        on the relative velocity of the cloud parcel and the fall-speed of the
        rain-droplets
        """
        rho_l = self.constants.rho_l

        # droplet-size distribution constant
        N0 = 1.0e7  # [m^-4]

        # size-distribtion length-scale
        l = (8.0 * rho_l * pi * N0 / (q_r * rho_c)) ** 0.25

        # fall-speed coefficient taken from the r > 0.5mm expression for
        # fall-speed from Herzog '98
        a_r = 201.0
        # reference density
        rho0 = 1.12

        # charateristic velocity, XXX: this is definitely wrong, but couldn't
        # work out how to do the integral at the time
        w_r = a_r * np.sqrt(1.0 / l * rho0 / rho_c)

        # rainout fraction
        f = w_r / w / self.l_pr

        return f * q_r

    def dFdz(self, z, F):
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

        # print "z={}m".format(z)

        # rho_e = self.environment.rho(z)
        T_e = self.environment.temp(z)
        try:
            qv_e = self.environment.rel_humidity(
                z
            ) * self.microphysics.parameterisations.pv_sat.qv_sat(T=T_e, p=p)
        except AttributeError:
            warnings.warn(
                "It seems the environmental profile doesn't define a relative humidity so we'll assume it's dry"
            )
            qv_e = 0.0
        qd_e = 1.0 - qv_e
        rho_e = self.cloud_mixture_density(
            p=p, T_c=T_e, qd_c=qd_e, qv_c=qv_e, ql_c=0.0, qr_c=0.0, qi_c=0.0
        )

        # calculate entrainment rate
        rho_c = self.cloud_mixture_density(
            p=p, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i, qr_c=q_r
        )
        g = self.constants.g
        B = (rho_e - rho_c) / rho_e
        mu = self._mu(r=r, w=w, B=B, z=z)

        # initiate gradient
        dFdz_ = np.zeros((Var.NUM,))

        # 1. Estimate change in vertical velocity with initial state
        dwdz_ = self.dw_dz(
            p=p,
            w_c=w,
            r_c=r,
            T_c=T,
            qd_c=q_d,
            qv_c=q_v,
            ql_c=q_l,
            qr_c=q_r,
            qi_c=q_i,
            rho_e=rho_e,
            mu=mu,
        )

        dFdz_[Var.w] = dwdz_

        # assume no condesates present in environment
        ql_e, qr_e, qi_e = 0.0, 0.0, 0.0

        # as well as tracer (water) changes from microphysics we have to consider entrainment
        dFdz_entrain__q = np.zeros((Var.NUM,))
        if self.entrain_hydrometeors:
            # dqd_dz__ent = mu/rho_c*(qd_e*rho_e - q_d*rho_c)
            dqd_dz__ent = mu * (qd_e - q_d)
            if dqd_dz__ent < 0.0:
                warnings.warn("WARNING: entrainment is moistening!")
            # dFdz_entrain__q[Var.q_v] = -q_v*dqd_dz__ent
            # dFdz_entrain__q[Var.q_l] = -q_l*dqd_dz__ent
            # dFdz_entrain__q[Var.q_r] = -q_r*dqd_dz__ent
            # dFdz_entrain__q[Var.q_i] = -q_i*dqd_dz__ent

            # this formulation gives a strong impact of entrainment on
            # hydrometeors, seems best I think
            dFdz_entrain__q[Var.q_v] = mu * (qv_e - q_v)
            dFdz_entrain__q[Var.q_l] = mu * (ql_e - q_l)
            dFdz_entrain__q[Var.q_r] = mu * (qr_e - q_r)
            dFdz_entrain__q[Var.q_i] = mu * (qi_e - q_i)

            # Old formulation: not exactly sure why this doesn't work, but the
            # above is from the perspective of entraining dry air instead of
            # detraining hydrometeors, maybe there's something important here?

            # dFdz_entrain__q[Var.q_v] = mu/rho_c*(qv_e*rho_e - q_v*rho_c)
            # dFdz_entrain__q[Var.q_l] = mu/rho_c*(ql_e*rho_e - q_l*rho_c)
            # dFdz_entrain__q[Var.q_r] = mu/rho_c*(qr_e*rho_e - q_r*rho_c)
            # dFdz_entrain__q[Var.q_i] = mu/rho_c*(qi_e*rho_e - q_i*rho_c)

            self.extra_vars.setdefault("dqv_dz__ent", []).append(
                dFdz_entrain__q[Var.q_v]
            )

        # 2. estimate new state from phase changes predicted by microphysics
        F[Var.p] = p  # make sure that pressure is set since the microphysics needs it
        dFdt_micro = self.microphysics.dFdt(
            F, t=z
        )  # passing in `z` as `t` doesn't make sense physically, but it will just be used for plotting later

        dFdz_micro = dFdt_micro / w  # w = dz/dt
        if self.entrain_liquid_static_energy or self.entrain_moist_static_energy:
            dFdz_micro[
                Var.T
            ] = 0.0  # temperature effect will be determined from temperature equation, microphysics simply provides gradients
        dFdz_ += dFdz_micro
        dFdz_ += dFdz_entrain__q

        # print "entrain:"
        # Var.print_formatted(dFdz_entrain__q)
        # print "micro:"
        # Var.print_formatted(dFdz_micro)

        # 3. Estimate temperature change forgetting about phase-changes for now (i.e. considering only adiabatic adjustment and entrainment)

        # in terms of the thermodynamics cloud-water and rain-water are treated identically
        dql_c__dz = dFdz_[Var.q_l] + dFdz_[Var.q_r]
        dqi_c__dz = dFdz_[Var.q_i]
        dqv_c__dz = dFdz_[Var.q_v]
        ql_c = q_l + q_r

        dTdz_s = self.dT_dz(
            r_c=r,
            T_c=T,
            qd_c=q_d,
            qv_c=q_v,
            ql_c=ql_c,
            qi_c=q_i,
            dql_c__dz=dql_c__dz,
            dqi_c__dz=dqi_c__dz,
            dqv_c__dz=dqv_c__dz,
            T_e=T_e,
            qv_e=qv_e,
            mu=mu,
        )

        dFdz_[Var.T] += dTdz_s
        dTdz_ = dFdz_[Var.T]

        # 4. Use post microphysics state (and phase changes from microphysics) to estimate radius change
        drdz_ = self.dr_dz(
            p=p,
            w_c=w,
            r_c=r,
            T_c=T,
            qd_c=q_d,
            qv_c=q_v,
            ql_c=q_l,
            qr_c=q_r,
            qi_c=q_i,
            dqv_c__dz=dqv_c__dz,
            dql_c__dz=dql_c__dz,
            dqi_c__dz=dqi_c__dz,
            dTc_dz=dTdz_,
            dw_dz=dwdz_,
            mu=mu,
        )

        dFdz_[Var.r] = drdz_

        # 5. Estimate fraction of rain that leaves parcel
        rho_c = self.cloud_mixture_density(
            p=p, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i, qr_c=q_r
        )
        dqr_dz__rainout = self.dqr_dz__rainout(rho_c=rho_c, q_r=q_r, w=w)
        dFdz_[Var.q_r] -= dqr_dz__rainout
        dFdz_[Var.q_pr] = dqr_dz__rainout

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
        model_desc = (
            r"FullSpecConcEqns ($C_D=%g$, $\beta=%g$, $l_{pr}=%s$), $\mu$-phys: %s"
            % (
                self.C_D,
                self.beta,
                ("\infty", "{:g}m".format(self.l_pr))[self.l_pr != np.inf],
                str(self.microphysics),
            )
        )

        if not self.temperature_dependent_latent_heats:
            model_desc += ", const $L_{heat}$"

        if not self.entrain_moist_static_energy:
            if self.entrain_liquid_static_energy:
                return model_desc + " LSE entr."
            else:
                return model_desc + " no therm. entr."
        else:
            return model_desc


class FullEquationsSatMicrophysics(FullThermodynamicsCloudEquations):
    """
    Same as the full model, but using the moist-adjustment microphysics is hardcoded in
    """

    def __init__(self, *args, **kwargs):
        if "microphysics" in kwargs:
            warnings.warn("%s is hardcoded to use the moist adjustment microphysics")
            kwargs["microphysics"] = cloud_microphysics.MoistAdjustmentMicrophysics()
        super(FullEquationsSatMicrophysics, self).__init__(*args, **kwargs)

    def dFdz(self, z, F):
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
        dwdz_ = self.dw_dz(
            p=p,
            w_c=w,
            r_c=r,
            T_c=T,
            qd_c=q_d,
            qv_c=q_v,
            ql_c=q_l,
            qi_c=q_i,
            rho_e=rho_e,
        )

        dFdz_[Var.w] = dwdz_

        # 2. Estimate temperature change forgetting about phase-changes for now (i.e. considering only adiabatic adjustment and entrainment)
        T_e = self.environment.temp(z)
        dTdz_s = self.dT_dz(
            r_c=r,
            T_c=T,
            qd_c=q_d,
            qv_c=q_v,
            ql_c=q_l,
            qi_c=q_i,
            dql_c__dz=0.0,
            dqi_c__dz=0.0,
            T_e=T_e,
        )

        F_s = np.copy(F)
        F_s[Var.T] += dTdz_s

        # 3. With new temperature estimate new state from phase changes predicted by microphysics
        F_s = self.microphysics._calc_adjusted_state(F_s, iterations=2)

        dTdz_ = F_s[Var.T] - F[Var.T]
        dFdz_[Var.T] = dTdz_

        dql_c__dz = F_s[Var.q_l] - F[Var.q_l]
        dqi_c__dz = F_s[Var.q_i] - F[Var.q_i]
        dFdz_[Var.q_v] = -dql_c__dz - dqi_c__dz
        dFdz_[Var.q_l] = dql_c__dz
        dFdz_[Var.q_i] = dqi_c__dz

        # 4. Use post microphysics state (and phase changes from microphysics) to estimate radius change
        drdz_ = self.dr_dz(
            p=p,
            w_c=w,
            r_c=r,
            T_c=T,
            qd_c=q_d,
            qv_c=q_v,
            ql_c=q_l,
            qi_c=q_i,
            dql_c__dz=dql_c__dz,
            dqi_c__dz=dqi_c__dz,
            dTc_dz=dTdz_,
            dw_dz=dwdz_,
        )

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

        rho_inv = (
            (qd_c * R_d + qv_c * R_v) * T_c / p + (ql_c + qr_c) / rho_l + qi_c / rho_i
        )

        return 1.0 / rho_inv

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

        R_s = (R_v * qv_c + R_d * qd_c) / (qd_c + qv_c)

        return p / (T_c * R_s)

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
        c_cm_p = cp_d * qd_c + cp_v * (qv_c + ql_c + qi_c)

        # heat capacity of environment mixture
        # qd_e, qv_e, ql_e, qi_e = self.qd_e, self.qv_e, self.ql_e, self.qi_e
        # TODO: Allow definition of moist environment
        qd_e, qv_e, ql_e, qi_e = 1.0, 0.0, 0.0, 0.0

        c_em_p = cp_d * qd_e + cp_v * (qv_e + ql_e + qi_e)

        # difference in environment and cloud moist static energy
        # XXX: There appears to be something wrong with the formulation that
        # includes liquid and ice, use just the liquid formulation for now
        # Ds = (c_em_p*T_e - ql_e*L_v - qi_e*L_s)\
        # -(c_cm_p*T_c - ql_c*L_v - qi_c*L_s)
        Ds = (c_em_p * T_e + ql_e * L_v) - (c_cm_p * T_c + ql_c * L_v)

        mu = self.beta / r_c

        return (
            -g / c_cm_p
            + mu * Ds / c_cm_p
            + L_v / c_cm_p * dql_c__dz
            + L_s / c_cm_p * dqi_c__dz
        )

    def dFdz(self, z, F):
        r = F[Var.r]
        w = F[Var.w]
        T = F[Var.T]
        q_v = F[Var.q_v]
        q_l = F[Var.q_l]
        q_i = F[Var.q_i]
        q_r = F[Var.q_r]
        q_d = 1.0 - q_v - q_l - q_i - q_r

        # cloud is assumed to be at same pressure as in environment
        T_e = self.environment.temp(z)

        dFdz_ = np.zeros((Var.NUM,))

        dFdz_[Var.w] = 0.0

        # 2. Estimate temperature change forgetting about phase-changes for now (i.e. considering only adiabatic adjustment and entrainment)
        dTdz_s = self.dT_dz(
            r_c=r,
            T_c=T,
            qd_c=q_d,
            qv_c=q_v,
            ql_c=q_l,
            qi_c=q_i,
            dql_c__dz=0.0,
            dqi_c__dz=0.0,
            T_e=T_e,
        )

        dFdz_[Var.T] = dTdz_s

        # 3. estimate new state from phase changes predicted by microphysics
        dFdt_micro = self.microphysics.dFdt(F, t=None)
        dFdz_micro = dFdt_micro / w  # w = dz/dt
        dFdz_ += dFdz_micro

        dFdz_[Var.r] = 0.0

        return dFdz_

    def __str__(self):
        return r"Fixed vertical velocity: w0={w0}m/s, $\mu$-phys: {mphys}".format(
            w0=self.w0, mphys=self.microphysics
        )
