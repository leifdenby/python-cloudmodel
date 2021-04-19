import numpy as np
import warnings


from ... import Var
from .. import methods as integration_methods
from ...plots import parcel as parcel_plots


class ParcelEvolution:
    def __init__(self, F, z, cloud_model, extra_vars={}):
        self.F = F
        self.z = z
        self.cloud_model = cloud_model
        self.extra_vars = extra_vars

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

        # solver = integration_methods.NewSolver(
        #     dFdz=self.cloud_model.dFdz, rtol=rtol, atol=atol, min_step=min_step
        # )

        # solver.set_initial_condition(initial_condition)
        # self.fail_on_nan = fail_on_nan

        # # TODO: work out how to extent state-vector to include microphysics variables
        # # mphys_extra_vars = {}
        # # self.microphysics.extra_vars = mphys_extra_vars

        # if stopping_criterion is None:
        #     stopping_criterion = self._stopping_criterion
        # F, z = solver.solve(z, stopping_criterion,)

        # if self._stopping_criterion(F_solution=F, z=z, k=len(z)-1):
        # F = F[:-1]
        # z = z[:-1]

        # for k, v in list(mphys_extra_vars.items()):
        # try:
        # mphys_extra_vars[k] = v[:-1]
        # except IndexError:
        # pass

        # mphys_extra_vars.update(self.extra_vars)
        rtol = 0.0
        solver = integration_methods.ScipyIntegrator(
            dFdz=self.cloud_model.dFdz, atol=atol, rtol=rtol
        )
        F, z = solver.solve(z=z, F0=initial_condition)

        return ParcelEvolution(F=F, z=z, cloud_model=self)

    def _stopping_criterion(self, z, F):
        if F[Var.T] > 300.0:
            print("Integration stopped: temperature got unphysically high")
            return True
        elif F[Var.w] < 0.0:
            print("Integration stopped: vertical velocity dropped to zero")
            print(F[Var.w])
            return True
        elif F[Var.r] > 10e3:
            print("Integration stopped: cloud radius became unreasonably high (r>10km)")
            return True
        elif z > 30e3:
            print("Integration stopped: height reached too high (z > 30km)")
            return True
        elif np.any(np.isnan(F)):
            print("Integration stopped: solution became nan")
            if self.fail_on_nan:
                raise Exception("Solution became nan")
            return True
        else:
            return False

    def _stopping_criterion_time(self, z, F):
        if hasattr(self, "stop_integration"):
            return True

        if F[Var.T] > 300.0:
            print("Integration stopped: temperature got unphysically high")
            return True
        elif F[Var.T] < 0.0:
            print("Integration stopped: temperature dropped below zero")
            return True
        elif F[Var.r] > 20e3:
            print("Integration stopped: cloud radius became unreasonably high (r>20km)")
            return True
        elif F[Var.z] < 0.0:
            print("Integration stopped: height below ground")
            return True
        elif F[Var.r] < 0.0:
            print("Integration stopped: radius dropped below zero")
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
            dFdt = dFdz * dzdt
            dFdt[Var.z] = dzdt

            # print "dFdz",
            # Var.print_formatted(dFdz)

            # print "dFdt=",
            # Var.print_formatted(dFdt)
            # print
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
