import numpy as np
import warnings


from ... import Var
from .. import methods as integration_methods
from ...plot import parcel as parcel_plots
from .stopping import DEFAULT_STOPPING_FUNCTIONS


class ParcelEvolution:
    def __init__(
        self, F, z, cloud_model, extra_vars={}, integration_stopping_reason=None
    ):
        self.F = F
        self.z = z
        self.cloud_model = cloud_model
        self.extra_vars = extra_vars
        self.integration_stopping_reason = integration_stopping_reason

    def plot(self, variables=("r", "w", "T")):
        return parcel_plots.plot_profiles(
            [
                self,
            ],
            variables=variables,
        )

    def __str__(self):
        return str(self.cloud_model)


class ParcelModelIntegrator(object):
    ALLOWS_FREEZING = False

    def __init__(self, cloud_model):
        self.cloud_model = cloud_model

    def __call__(
        self,
        initial_condition,
        z,
        stopping_criterion=None,
        tolerance=1e-3,
        fail_on_nan=False,
        method="RK45",
    ):
        self._validate_initial_state(initial_condition)

        atol_q = 1.0e-3 * 1.0e-4
        T_tol = 0.01
        atol = Var.make_state(
            p=100.0,
            T=T_tol,
            q_v=atol_q,
            q_l=atol_q,
            w=0.1,
            r=0.1,
            q_r=atol_q,
            z=10.0,
            q_i=atol_q,
            q_pr=atol_q,
        )

        rtol = 0.0
        solver = integration_methods.ScipyIntegrator(
            dFdz=self.cloud_model.dFdz,
            atol=atol,
            rtol=rtol,
            stopping_functions=DEFAULT_STOPPING_FUNCTIONS,
            method=method,
        )
        F, z, stopping_reason = solver.solve(z=z, F0=initial_condition)

        return ParcelEvolution(
            F=F,
            z=z,
            cloud_model=self.cloud_model,
            integration_stopping_reason=stopping_reason,
        )

    def dFdt(self, F, t):
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

    @staticmethod
    def _validate_initial_state(F):
        if F[Var.T] == 0.0:
            raise Exception("Should have temperature > 0.0")
        if F[Var.p] == 0.0:
            raise Exception("Should have pressure > 0.0")

    def integrate_in_time(
        self, initial_condition, t, SolverClass, stopping_criterion=None, tolerance=1e-3
    ):
        warnings.warn("Integration in time is still experimental")

        self._validate_initial_state(initial_condition)

        solver = SolverClass(
            self.dFdt,
            rtol=0.0,
            atol=tolerance,
        )
        solver.set_initial_condition(np.array(initial_condition))

        if stopping_criterion is None:
            stopping_criterion = self._stopping_criterion_time

        F, t = solver.solve(t, stopping_criterion)

        return ParcelEvolution(F=F, z=F[:, Var.z], cloud_model=self)
