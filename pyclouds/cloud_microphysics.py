"""
Collection of microphysics routines for use with cloud-model integration.
"""
import warnings
import numpy as np
import odespy
from scipy.constants import pi

from pyclouds.cloud_equations import Var
from common import AttrDict
from pyclouds.utils import default_constants
from pyclouds.plotting import plot_hydrometeor_evolution
from pyclouds import parameterisations as default_parameterisations

# try:
    # from ccfm.ccfmfortran import microphysics as ccfm_microphysics
# except ImportError:
    # # import pure python version instead
from ccfm.ccfmpython import microphysics as ccfm_microphysics

class HydrometeorEvolution:
    def __init__(self, F, t, model, integration_kwargs={}):
        self.F = F
        self.t = t
        self.model = model
        self.integration_kwargs = integration_kwargs

    def plot(self):
        plot_hydrometeor_evolution([self,])

    def __str__(self):
        s = str(self.model)
        if len(self.integration_kwargs) > 0:
            s += " (%s)" % ", ".join(["%s: %s" % (k, str(v)) for (k, v) in self.integration_kwargs.items()])
        return s


class BaseMicrophysicsModel(object):
    def __init__(self, constants=default_constants, *args, **kwargs):
        self.parameterisations = kwargs.get('parameterisations', default_parameterisations)
        self.constants = AttrDict(constants)
        if not hasattr(constants, 'R_d'):
            self.constants.R_d = self.constants.cp_d - self.constants.cv_d
            self.constants.R_v = self.constants.cp_v - self.constants.cv_v

    def dFdt(F, t):
        raise NotImplemented

    def _stopping_criterion(self, F_solution, t, k):
        return False

    def integrate(self, initial_condition, t, p0, SolverClass=odespy.RKFehlberg, stopping_criterion=None, tolerance=1e-3):
        # When integrating the microphysics on it's own we assume constant
        # pressure which must be provided when integrating
        self.p0 = p0
        derriv_f = lambda F, t: self.dFdt(F=F, p=p0, t=t)
        solver = SolverClass(derriv_f, rtol=0.0, atol=tolerance,)
        solver.set_initial_condition(initial_condition)

        if stopping_criterion is None:
            stopping_criterion=self._stopping_criterion
        F, t = solver.solve(t, stopping_criterion,)

        return HydrometeorEvolution(F=F, t=t, model=self)

class DummyMicrophysics(BaseMicrophysicsModel):
    """
    Dummy microphysics implementation that doesn't affect the state at all.
    """


    def dFdt(self, F, p):
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
    def __init__(self, *args, **kwargs):
        super(FiniteCondensationTimeMicrophysics, self).__init__(*args, **kwargs)

        self.N0 = 200*1.e6  # initial aerosol number concentration [m-3]
        self.r0 = 0.1e-6  # cloud droplet initial radius
        self.r_crit = 5.0e-6  # critical cloud droplet radius [m] after which the number of cloud droplets is increased

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
        rho_l = self.constants.rho_l
        Lv = self.constants.L_v
        R_v = self.constants.R_v
        R_d = self.constants.R_d
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


        r_c = (ql*rho/(4./3.*pi*self.N0*rho_l))

        if r_c > self.r_crit:
            # if cloud droplet radius with initial number of droplets is larger
            # than a critial size assume that instead more droplets are made,
            # all with the critical radius
            Nc = q_c*rho/(4./3.*pi*rho_l*self.r_crit**3.0)
            r_c = self.r_crit
        else:
            Nc = self.N0
            r_c = self.r0


        # condensation evaporation of cloud droplets (given number of droplets
        # and droplet radius calculated above)

        # TODO: This needs double checking
        Sw = qv/self.qv_sat(T=T, p=p)

        # TODO: move the thermal conductivity of and diffusion through air somewhere else
        Ka = self.parameterisations.Ka(T=T)
        Fk = (Lv/(R_v*T) - 1)*Lv/(Ka*T)

        pv_sat = self.parameterisations.pv_sat(T=T)
        Dv = self.parameterisations.Dv(T=T)
        Fd = R_v*T/(pv_sat*Dv)

        # compute rate of change of condensate from diffusion
        dqc_dt = 4*pi*1./rho*Nc*r_c*(Sw - 1.0)/(Fk + Fd)

        dFdz = np.zeros((Var.NUM))
        dFdz[Var.q_l] = dqc_dt
        dFdz[Var.q_v] = -dqc_dt
        dFdz[Var.T] = Lv/cp_m*dqc_dt

        # print dqc_dt, Sw
        # import ipdb
        # ipdb.set_trace()

        return dFdz


    def __call__(self, F, p):
        pass

    def __str__(self):
        return "Finite condensation time"
