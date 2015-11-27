"""
Collection of microphysics routines for use with cloud-model integration.
"""
import warnings
import numpy as np
import odespy
from scipy.constants import pi

from pyclouds.cloud_equations import Var
from common import AttrDict, default_constants, make_related_constants, ATHAM_constants
from pyclouds.plotting import plot_hydrometeor_evolution
from pyclouds import parameterisations

try:
    import unified_microphysics.fortran as unified_microphysics
    from unified_microphysics.tests.test_common import um_constants
except ImportError:
    unified_microphysics = None

# try:
    # from ccfm.ccfmfortran import microphysics as ccfm_microphysics
# except ImportError:
    # # import pure python version instead
from ccfm.ccfmpython import microphysics as ccfm_microphysics

class PyCloudsUnifiedMicrophysicsStateMapping():
    """
    Utility for providing mapping between pycloud state representation and
    unified-microphysics model representation.

    TODO: Make pyclouds use indexing from `unified microphysics`
    
    XXX: this is kinda hacky, I only reassign state variables that are actually
    used in the initialized model in `unified microphysics`
    """
    def __init__(self):
        register = unified_microphysics.microphysics_register
        # Fortran indexing starts at 1
        self.idx_water_vapour = register.idx_water_vapour-1
        self.idx_cwater = register.idx_cwater-1
        self.idx_rain = register.idx_rain-1
        self.idx_cice = register.idx_cice-1
        self.idx_temp = register.idx_temp-1
        self.idx_pressure = register.idx_pressure-1

        self.n_vars = register.n_variables

    def um_pycloud(self, y):
        F = np.zeros((Var.NUM))
        if self.idx_water_vapour != -1:
            F[Var.q_v] = y[self.idx_water_vapour]
        if self.idx_cwater != -1:
            F[Var.q_l] = y[self.idx_cwater]
        if self.idx_rain != -1:
            F[Var.q_r] = y[self.idx_rain]
        if self.idx_cice != -1:
            F[Var.q_i] = y[self.idx_cice]
        F[Var.T] = y[self.idx_temp]
        return F, y[self.idx_pressure]

    def pycloud_um(self, F, p):
        y = np.zeros((self.n_vars,))

        if self.idx_water_vapour != -1:
            y[self.idx_water_vapour] = F[Var.q_v]
        if self.idx_cwater != -1:
            y[self.idx_cwater] = F[Var.q_l]
        if self.idx_rain != -1:
            y[self.idx_rain] = F[Var.q_r]
        if self.idx_cice != -1:
            y[self.idx_cice] = F[Var.q_i]
        y[self.idx_temp] = F[Var.T]
        y[self.idx_pressure] = p
        return y

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
        F_adjusted = self._calc_adjusted_state(F=F0, p=p0, iterations=iterations)

        if F_adjusted[Var.q_l] < 0.0:
            F[1:] = F0
        else:
            F[1:] = F_adjusted

        return HydrometeorEvolution(F=F, t=t, model=self, integration_kwargs={ 'iterations': iterations, })

    def __call__(self, *args, **kwargs):
        raise Exception("Since moist adjustment is instaneous it cannot be called to provided a differential")

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
        self.disable_rain = kwargs.get('disable_rain', False)

    def calc_mixture_density(self, qd, qv, ql, qi, qr, p, T):
        warnings.warn("EoS calculation stored within microphysics, should really use something defined externally")

        R_d = self.constants.R_d
        R_v = self.constants.R_v
        rho_l = self.constants.rho_l
        rho_i = self.constants.rho_i

        rho_inv = (qd*R_d + qv*R_v)*T/p + (ql+qr)/rho_l + qi/rho_i
        
        return 1.0/rho_inv

    def qv_sat(self, T, p):
        return self.parameterisations.pv_sat.qv_sat(T=T, p=p)

    def cp_m(self, F):
        qv = F[Var.q_v]
        ql = F[Var.q_l]
        qr = F[Var.q_r]
        qi = F[Var.q_i]
        qd = 1. - qv - ql - qr - qi

        cp_d = self.constants.cp_d
        cp_v = self.constants.cp_v
        cp_l = self.constants.cp_l

        if qi > 0.0:
            raise NotImplementedError

        return cp_d*qd + (ql + qr)*cp_l + qv*cp_v

    def cv_m(self, F):
        qv = F[Var.q_v]
        ql = F[Var.q_l]
        qr = F[Var.q_r]
        qi = F[Var.q_i]
        qd = 1. - qv - ql - qr - qi

        cv_d = self.constants.cv_d
        cv_v = self.constants.cv_v
        cv_l = self.constants.cv_l

        if qi > 0.0:
            raise NotImplementedError

        return cv_d*qd + (ql + qr)*cv_l + qv*cv_v

    def dFdt(self, F, t, p):
        Lv = self.constants.L_v
        cp_d = self.constants.cp_d
        cp_v = self.constants.cp_v

        qv = F[Var.q_v]
        ql = F[Var.q_l]
        qr = F[Var.q_r]
        qi = F[Var.q_i]

        # the integrator may have put us below zero, this is non-conservative,
        # but I don't have any other options here since I can't modify the
        # integrator's logic
        if ql < 0.0:
            ql = 0.0
            F[Var.q_l] = ql

        qd = 1. - qv - ql - qr - qi
        T = F[Var.T]

        cp_m = self.cp_m(F=F)

        # mixture density
        rho = self.calc_mixture_density(qd=qd, qv=qv, ql=ql, qi=qi, qr=qr, p=p, T=T)
        # gas density
        rho_g = self.calc_mixture_density(qd=qd, qv=qv, ql=0., qi=0., qr=0., p=p, T=T)

        dql_dt = self.dql_dt__cond_evap(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=T, p=p)

        qg = qv + qd
        dqr_dt_1 = self._dqr_dt__autoconversion(ql=ql, qg=qg, rho_g=rho_g)
        dqr_dt_2 = self._dqr_dt__accretion(ql=ql, qg=qg, qr=qr, rho_g=rho_g)

        dqr_dt = dqr_dt_1 + dqr_dt_2

        if self.disable_rain:
            dqr_dt = 0.0

        dFdz = np.zeros((Var.NUM))
        dFdz[Var.q_l] =  dql_dt -dqr_dt
        dFdz[Var.q_v] = -dql_dt
        dFdz[Var.q_r] =          dqr_dt
        dFdz[Var.T] = Lv/cp_m*dql_dt

        if self.debug and hasattr(self, 'extra_vars'):
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

        return max(dqr_dt, 0.0)


    def dql_dt__cond_evap(self, rho, rho_g, qv, ql, T, p):
        Lv = self.constants.L_v
        R_v = self.constants.R_v
        R_d = self.constants.R_d
        rho_l = self.constants.rho_l

        # condensation evaporation of cloud droplets (given number of droplets
        # and droplet radius calculated above)
        Sw = qv/self.qv_sat(T=T, p=p)

        if ql == 0.0 and Sw > 1.0:
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


        Ka = self.parameterisations.Ka(T=T)
        Fk = (Lv/(R_v*T) - 1)*Lv/(Ka*T)

        pv_sat = self.parameterisations.pv_sat(T=T)
        Dv = self.parameterisations.Dv(T=T, p=p)
        Fd = R_v*T/(pv_sat*Dv)

        # compute rate of change of condensate from diffusion
        dql_dt = 4*pi*1./rho*Nc*r_c*(Sw - 1.0)/(Fk + Fd)

        if self.debug and hasattr(self, 'extra_vars'):
            self.extra_vars.setdefault('r_c', []).append(r_c)
            self.extra_vars.setdefault('Nc', []).append(Nc)

        return dql_dt


    def __str__(self):
        s = "Finite condensation rate ($r_{crit}=%gm$)" % (self.r_crit)
        if self.disable_rain:
            s += " without rain"
        return s

class FC_min_radius(FiniteCondensationTimeMicrophysics):

    def dql_dt__cond_evap(self, rho, rho_g, qv, ql, T, p):
        Lv = self.constants.L_v
        R_v = self.constants.R_v
        R_d = self.constants.R_d
        rho_l = self.constants.rho_l

        # condensation evaporation of cloud droplets (given number of droplets
        # and droplet radius calculated above)
        Sw = qv/self.qv_sat(T=T, p=p)

        r_c = (ql*rho/(4./3.*pi*self.N0*rho_l))**(1./3.)

        if r_c > self.r_crit:
            # if cloud droplet radius with initial number of droplets is larger
            # than a critial size assume that instead more droplets are made,
            # all with the critical radius
            Nc = ql*rho/(4./3.*pi*rho_l*self.r_crit**3.0)
            r_c = self.r_crit
        elif r_c < self.r0:
            r_c = self.r0
            Nc = self.N0
        else:
            Nc = self.N0

        Ka = self.parameterisations.Ka(T=T)
        Fk = (Lv/(R_v*T) - 1)*Lv/(Ka*T)

        pv_sat = self.parameterisations.pv_sat(T=T)
        Dv = self.parameterisations.Dv(T=T, p=p)
        Fd = R_v*T/(pv_sat*Dv)

        # compute rate of change of condensate from diffusion
        dql_dt = 4*pi*1./rho*Nc*r_c*(Sw - 1.0)/(Fk + Fd)

        if self.debug and hasattr(self, 'extra_vars'):
            self.extra_vars.setdefault('r_c', []).append(r_c)
            self.extra_vars.setdefault('Nc', []).append(Nc)

        return dql_dt

    def __str__(self):
        return super(FC_min_radius, self).__str__() + r", $r_{min}=%gm$" % self.r0

class FortranNoIceMicrophysics(BaseMicrophysicsModel):
    """
    Wrapper for state variable differentials calculated in
    `unified-microphysics` Fortran library.
    """

    def __init__(self, *args, **kwargs):
        super(FortranNoIceMicrophysics, self).__init__(*args, **kwargs)
        if unified_microphysics is None:
            raise Exception("Couldn't import the `unified_microphysics` library, please symlink it to `unified_microphysics.so`")

        if not 'model_constraint' in kwargs:
            warnings.warn("`model_constraint` not provided, assuming isometric")
        model_constraint = kwargs.get('model_constraint', 'isometric')

        self.model_constraint = model_constraint

        unified_microphysics.microphysics_pylib.init('no_ice', model_constraint)

        constants = um_constants

        self.parameterisations = parameterisations.ParametersationsWithSpecificConstants(constants=constants)
        self.qv_sat = self.parameterisations.pv_sat.qv_sat

    def dFdt(self, F, t, p):
        state_mapping = PyCloudsUnifiedMicrophysicsStateMapping()

        ql = F[Var.q_l]
        # the integrator may have put us below zero, this is non-conservative,
        # but I don't have any other options here since I can't modify the
        # integrator's logic
        if ql < 0.0:
            ql = 0.0
            F[Var.q_l] = ql

        y = state_mapping.pycloud_um(F, p)

        c_m = unified_microphysics.microphysics_pylib.mixture_heat_capacity(y)

        dydt = unified_microphysics.mphys_no_ice.dydt(y=y, t=0.0, c_m=c_m)

        dFdz, _ = state_mapping.um_pycloud(y=dydt)

        return dFdz

    def __str__(self):
        return "Fortran 'no_ice' model (%s) with odespy integrator" % self.model_constraint

class ExplicitFortranModel:
    """
    Calls the integrator implemented in Fortran in um
    """

    def __init__(self, *args, **kwargs):
        if unified_microphysics is None:
            raise Exception("Couldn't import the `unified_microphysics` library, please symlink it to `unified_microphysics.so`")

        constants = um_constants

        unified_microphysics.microphysics_pylib.init('no_ice')

        self.parameterisations = parameterisations.ParametersationsWithSpecificConstants(constants=constants)
        self.qv_sat = self.parameterisations.pv_sat.qv_sat

    def integrate(self, initial_condition, t, p0, SolverClass=odespy.RKFehlberg, stopping_criterion=None, tolerance=1e-3):
        # When integrating the microphysics on it's own we assume constant
        # pressure which must be provided when integrating
        self.p0 = p0

        state_mapping = PyCloudsUnifiedMicrophysicsStateMapping()

        y = state_mapping.pycloud_um(initial_condition, p0)
        F = [state_mapping.um_pycloud(y=y)[0],]
        t_ = [t[0],]

        for tn in range(len(t)-1):
            # modifies `y` in-place
            try:
                unified_microphysics.microphysics_integration.integrate_with_message(y=y, t0=t[tn], t_end=t[tn+1])
                F.append(state_mapping.um_pycloud(y=y)[0])
                t_.append(t[tn+1])
            except Exception as e:
                print "(%s): Integration stopped, %s" % (str(self), str(e))
                break

        return HydrometeorEvolution(F=np.array(F), t=np.array(t_), model=self,)


    def __str__(self):
        return "pure fortran `no_ice` model with rkf34"

class OldATHAMKesslerFortran:
    """
    Uses refactored implementation of "Kessler microphysics" from ATHAM
    implemented in the `unified microphysics` codebase.
    """

    def __init__(self, *args, **kwargs):
        if unified_microphysics is None:
            raise Exception("Couldn't import the `unified_microphysics` library, please symlink it to `unified_microphysics.so`")

        unified_microphysics.microphysics_pylib.init('kessler_old', 'isobaric')

        self.parameterisations = parameterisations.ParametersationsWithSpecificConstants(constants=ATHAM_constants)
        self.qv_sat = self.parameterisations.pv_sat.qv_sat

    def integrate(self, initial_condition, t, p0):
        self.p0 = p0

        state_mapping = PyCloudsUnifiedMicrophysicsStateMapping()

        y = state_mapping.pycloud_um(initial_condition, p0)
        F = [state_mapping.um_pycloud(y=y)[0],]
        t_ = [t[0],]

        for tn in range(len(t)-1):
            # modifies `y` in-place
            try:
                unified_microphysics.mphys_kessler_old.integrate_isobaric(y=y, t0=t[tn], t_end=t[tn+1])
                F_new = state_mapping.um_pycloud(y=y)[0]
                F.append(F_new)
                t_.append(t[tn+1])
                if F_new[Var.T] < 0.0:
                    raise ValueError("Temperature became negative")
            except Exception as e:
                print "(%s): Integration stopped, %s" % (str(self), str(e))
                break


        return HydrometeorEvolution(F=np.array(F), t=np.array(t_), model=self,)


    def __str__(self):
        return 'old ATHAM "Kessler microphysics"'
