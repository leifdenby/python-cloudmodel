"""
Collection of microphysics routines for use with cloud-model integration.
"""
import warnings
import numpy as np
import odespy
from scipy.constants import pi

from pyclouds.cloud_equations import Var
from common import AttrDict, default_constants, make_related_constants
from pyclouds.plotting import plot_hydrometeor_evolution
from pyclouds import parameterisations

try:
    import unified_microphysics
except ImportError:
    unified_microphysics = None

# try:
    # from ccfm.ccfmfortran import microphysics as ccfm_microphysics
# except ImportError:
    # # import pure python version instead
from ccfm.ccfmpython import microphysics as ccfm_microphysics

class HydrometeorEvolution:
    def __init__(self, F, t, model, integration_kwargs={}, extra_vars={}):
        self.F = F
        self.t = t
        self.model = model
        self.integration_kwargs = integration_kwargs
        self.extra_vars = extra_vars

    def plot(self):
        plot_hydrometeor_evolution([self,])

    def __str__(self):
        s = str(self.model)
        if len(self.integration_kwargs) > 0:
            s += " (%s)" % ", ".join(["%s: %s" % (k, str(v)) for (k, v) in self.integration_kwargs.items()])
        return s

class BaseMicrophysicsModel(object):
    def __init__(self, constants=default_constants, *args, **kwargs):
        constants = make_related_constants(constants)

        self.parameterisations = parameterisations.ParametersationsWithSpecificConstants(constants=constants)
        self.constants = AttrDict(constants)

    def dFdt(F, t):
        raise NotImplemented

    def _stopping_criterion(self, F_solution, t, k):
        return False

    def integrate(self, initial_condition, t, p0, SolverClass=odespy.RKFehlberg, stopping_criterion=None, tolerance=1e-3):
        # When integrating the microphysics on it's own we assume constant
        # pressure which must be provided when integrating
        self.p0 = p0
        self.extra_vars = {}
        derriv_f = lambda F, t: self.dFdt(F=F, p=p0, t=t)
        solver = SolverClass(derriv_f, rtol=0.0, atol=tolerance,)
        solver.set_initial_condition(initial_condition)

        if stopping_criterion is None:
            stopping_criterion=self._stopping_criterion
        F, t = solver.solve(t, stopping_criterion,)

        return HydrometeorEvolution(F=F, t=t, model=self, extra_vars=self.extra_vars,)

    def __call__(self, F, p, dt):
        return self.dFdt(F=F, p=p)*dt

class DummyMicrophysics(BaseMicrophysicsModel):
    """
    Dummy microphysics implementation that doesn't affect the state at all.
    """

    def dFdt(self, F, p):
        return np.zeros((Var.NUM))

    def __call__(self, F, p, dt):
        return F

class MoistAdjustmentMicrophysics(BaseMicrophysicsModel):
    """
    Adjust temperature and specific concentration of water vapour to saturation
    point, assuming pressure constant.

    Uses `moist_adjust` method from `mo_ccfm_cloudbase.f90`
    """

    def integrate(self, initial_condition, t, p0, iterations=1, *args, **kwargs):
        """
        Fake integrate method, moist adjustment instantenously turns all supersaturation into condensed water.
        """

        self.p0 = p0

        F0 = initial_condition
        F = np.zeros((len(t), Var.NUM))
        F[0] = initial_condition
        F[1:] = self._calc_adjusted_state(F=F0, p=p0, iterations=iterations)

        return HydrometeorEvolution(F=F, t=t, model=self, integration_kwargs={ 'iterations': iterations, })

    def __call__(self, F, p, dt):
        return self._calc_adjusted_state(F=F, p=p, iterations=3)

    def qv_sat(self, T, p):
        return self.parameterisations.pv_sat.qv_sat(T=T, p=p)

    def _calc_adjusted_state(self, F, p, iterations):
        """
        Calculate approximation saturation state using first-order Taylor
        expansion of the moist enthalpy equation at constant pressure:

        cp*dT = -L_v*dq

        Following suggestions in notes by Peter Bechtold
        """
        cp_d = self.constants.cp_d
        cp_v = self.constants.cp_v
        L_v = self.constants.L_v

        qv = F[Var.q_v]
        ql = F[Var.q_l]
        qr = F[Var.q_r]
        qi = F[Var.q_i]
        qd = 1. - qv - ql - qr - qi

        T = F[Var.T]
        qv = F[Var.q_v]

        for n in range(iterations):
            cp_m = cp_d*qd + (ql + qr + qi + qv)*cp_v
            qv_sat = self.parameterisations.pv_sat.qv_sat(T=T, p=p)
            dqv_sat__dT = self.parameterisations.pv_sat.dqv_sat__dT(T=T, p=p)

            dT = L_v/cp_m*(qv - qv_sat)/(1 + L_v/cp_m*dqv_sat__dT)

            T = T+dT
            qv = qv - cp_m/L_v*dT
            ql = ql + cp_m/L_v*dT

        Fs = np.copy(F)
        Fs[Var.q_v] = qv
        Fs[Var.q_l] = ql
        Fs[Var.T] = T

        return Fs

    def __str__(self):
        return "Moist adjustment"

class FiniteCondensationTimeMicrophysics(BaseMicrophysicsModel):
    def __init__(self, r_crit=5.0e-6, *args, **kwargs):
        super(FiniteCondensationTimeMicrophysics, self).__init__(*args, **kwargs)

        self.N0 = 200*1.e6  # initial aerosol number concentration [m-3]
        self.r0 = 0.1e-6  # cloud droplet initial radius
        self.r_crit = r_crit  # critical cloud droplet radius [m] after which the number of cloud droplets is increased
        self.debug = True

    def _calc_mixture_density(self, qd, qv, ql, qi, qr, p, T):
        warnings.warn("EoS calculation stored within microphysics, should really use something defined externally")

        R_d = self.constants.R_d
        R_v = self.constants.R_v
        rho_l = self.constants.rho_l
        rho_i = self.constants.rho_i

        rho_inv = (qd*R_d + qv*R_v)*T/p + (ql+qr)/rho_l + qi/rho_i
        
        return 1.0/rho_inv

    def qv_sat(self, T, p):
        return self.parameterisations.pv_sat.qv_sat(T=T, p=p)

    def dFdt(self, F, t, p):
        Lv = self.constants.L_v
        cp_d = self.constants.cp_d
        cp_v = self.constants.cp_v

        qv = F[Var.q_v]
        ql = F[Var.q_l]
        qr = F[Var.q_r]
        qi = F[Var.q_i]
        qd = 1. - qv - ql - qr - qi
        T = F[Var.T]

        cp_m = cp_d*qd + (ql + qr + qi + qv)*cp_v

        # mixture density
        rho = self._calc_mixture_density(qd=qd, qv=qv, ql=ql, qi=qi, qr=qr, p=p, T=T)
        # gas density
        rho_g = self._calc_mixture_density(qd=qd, qv=qv, ql=0., qi=0., qr=0., p=p, T=T)

        dql_dt = self._dql_dt__cond_evap(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=T, p=p)

        qg = qv + qd
        dqr_dt_1 = self._dqr_dt__autoconversion(ql=ql, qg=qg, rho_g=rho_g)
        dqr_dt_2 = self._dqr_dt__accretion(ql=ql, qg=qg, qr=qr, rho_g=rho_g)

        dqr_dt = dqr_dt_1 + dqr_dt_2

        dFdz = np.zeros((Var.NUM))
        dFdz[Var.q_l] =  dql_dt -dqr_dt
        dFdz[Var.q_v] = -dql_dt
        dFdz[Var.q_r] =          dqr_dt
        dFdz[Var.T] = Lv/cp_m*dql_dt

        if self.debug:
            self.extra_vars.setdefault('t_substeps', []).append(t)

        return dFdz

    def _dqr_dt__autoconversion(self, ql, qg, rho_g):
        """
        Create rain droplets through collision and coalescence of cloud
        droplets
        """
        k_c = 1.e-3
        a_c = 5.e-4

        dqr_dt = k_c*(ql - qg/rho_g*a_c)

        # TODO: Is it reasonable to limit this? It seems the autoconversion
        # equation doesn't really make sense to describe breakup
        dqr_dt = max(0., dqr_dt)

        return dqr_dt

    def _dqr_dt__accretion(self, rho_g, ql, qg, qr):
        rho_l = self.constants.rho_l
        G3p5 = 3.32399614155  # = Gamma(3.5)
        N0r = 1.e7  # [m^-4]
        a_r = 201.0  # [m^.5 s^-1]
        rho0 = 1.12

        lambda_r = (pi*(qg*rho_l)/(qr*rho_g)*N0r)**(1./4.)

        dqr_dt = pi/4.*N0r*a_r*np.sqrt(rho0/rho_g)*G3p5*lambda_r**(-3.5)*ql

        return 0.*max(dqr_dt, 0.0)


    def _dql_dt__cond_evap(self, rho, rho_g, qv, ql, T, p):
        Lv = self.constants.L_v
        R_v = self.constants.R_v
        R_d = self.constants.R_d
        rho_l = self.constants.rho_l

        if ql == 0.0:
            r_c = self.r0
        else:
            r_c = (ql*rho/(4./3.*pi*self.N0*rho_l))**(1./3.)

        if r_c > self.r_crit:
            # if cloud droplet radius with initial number of droplets is larger
            # than a critial size assume that instead more droplets are made,
            # all with the critical radius
            Nc = ql*rho/(4./3.*pi*rho_l*self.r_crit**3.0)
            r_c = self.r_crit
        else:
            Nc = self.N0

        # condensation evaporation of cloud droplets (given number of droplets
        # and droplet radius calculated above)
        Sw = qv/self.qv_sat(T=T, p=p)

        Ka = self.parameterisations.Ka(T=T)
        Fk = (Lv/(R_v*T) - 1)*Lv/(Ka*T)

        pv_sat = self.parameterisations.pv_sat(T=T)
        Dv = self.parameterisations.Dv(T=T)
        Fd = R_v*T/(pv_sat*Dv)

        # compute rate of change of condensate from diffusion
        dql_dt = 4*pi*1./rho*Nc*r_c*(Sw - 1.0)/(Fk + Fd)

        if self.debug:
            self.extra_vars.setdefault('r_c', []).append(r_c)
            self.extra_vars.setdefault('Nc', []).append(Nc)

        return dql_dt


    def __call__(self, F, p):
        pass

    def __str__(self):
        return "Finite condensation time ($r_{crit}=%gm$)" % (self.r_crit)

class FortranNoIceMicrophysics(BaseMicrophysicsModel):
    """
    Wrapper for state variable differentials calculated in
    `unified-microphysics` Fortran library.
    """

    def __init__(self, *args, **kwargs):
        super(FortranNoIceMicrophysics, self).__init__(*args, **kwargs)
        if unified_microphysics is None:
            raise Exception("Couldn't import the `unified_microphysics` library, please symlink it to `unified_microphysics.so`")

        unified_microphysics.microphysics_pylib.init('no_ice')

    def dFdt(self, F, t, p):
        # TODO: construct state variable structure using information from
        # fortran library, instead of mapping explicitly here

        register = unified_microphysics.microphysics_register
        # Fortran indexing starts at 1
        idx_water_vapour = register.idx_water_vapour-1
        idx_cwater = register.idx_cwater-1
        idx_rain = register.idx_rain-1
        idx_cice = register.idx_cice-1

        n_gases = register.n_gases
        n_moments__max = register.n_moments__max
        n_solids = register.n_solids

        q_g = np.zeros((n_gases,))
        q_tr = np.zeros((n_solids, n_moments__max))

        q_g[idx_water_vapour] = F[Var.q_v]
        q_tr[idx_cwater,0] = F[Var.q_l]
        q_tr[idx_rain,0] = F[Var.q_r]
        q_tr[idx_cice,0] = F[Var.q_i]
        T = F[Var.T]

        dqdt_g, dqdt_tr, _ = unified_microphysics.microphysics_pylib.dqdt(q_g=q_g, q_tr=q_tr, temp=T, pressure=p)

        dFdz = np.zeros((Var.NUM))
        dFdz[Var.q_v] = dqdt_g[idx_water_vapour]
        dFdz[Var.q_l] = dqdt_tr[idx_cwater]
        dFdz[Var.q_r] = dqdt_tr[idx_rain]
        dFdz[Var.q_i] = dqdt_tr[idx_cice]

        return dFdz

    def __str__(self):
        return "Fortran ('no_ice') model"
