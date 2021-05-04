import sys
import numpy as np
from scipy.integrate import solve_ivp

from .. import Var


class ScipyIntegrator:
    """
    Wrapper for calling scipy's ODE integration so that we can pass in stopping
    criteria and check the output easily
    """

    def __init__(self, dFdz, atol, rtol, method="RK45", stopping_functions=[]):
        self.dFdz = dFdz
        self.atol = atol
        self.rtol = rtol
        self.debug = False
        [setattr(fn, "terminal", True) for fn in stopping_functions]
        self.stopping_functions = stopping_functions
        self.method = method

    def solve(self, z, F0):
        res = solve_ivp(
            fun=self.dFdz,
            method=self.method,
            t_span=[z.min(), z.max()],
            y0=F0,
            atol=self.atol,
            events=self.stopping_functions,
        )

        F, z = res.y.T, res.t

        if res.status == 0:
            return F, z, None
        elif res.status == 1:
            # a stopping criterion stopped the integration, find out which one(s)
            function_names = [fn.__name__ for fn in self.stopping_functions]
            triggered_fns = {}
            for i, z_triggers in enumerate(res.t_events):
                if len(z_triggers) == 0:
                    pass
                elif len(z_triggers) == 1:
                    triggered_fns[function_names[i]] = z_triggers[0]
                else:
                    raise NotImplementedError(len(z_triggers))

            stopping_reason = ", ".join(
                [f"{k} triggered at {v}" for (k, v) in triggered_fns.items()]
            )

            return F, z, stopping_reason
        else:
            stopping_reason = res.message
            return F, z, stopping_reason


class NewSolver:
    """
    Simple Runge-Kutta 4-5 predictor-corrector integrator
    """

    def __init__(self, dFdz, atol, rtol, stopping_func):
        self.dFdz = dFdz
        self.atol = atol
        self.rtol = rtol
        self.min_step = 0.01
        self.debug = False
        self.min_step = 0.01
        self.stopping_func = stopping_func

    def solve(self, z, F0):
        f = self.dFdz
        x0 = F0
        a = z.min()
        b = z.max()
        t = a
        x = x0
        h = 1.0

        # Initialize arrays that will be returned

        T = np.array([t])
        X = np.array([x])

        S = np.array(
            [
                np.nan,
            ]
        )

        kk = 0

        while t < b:
            kk += 1
            # Adjust step size when we get to last interval

            if self.stopping_func(X, T, len(T) - 1):
                return (X, T)

            if self.debug:
                print()
                print("kk, t, h", kk, t, h)

                print("start")
                Var.print_formatted(x)
            else:
                if int(t) % 100 == 0:
                    print("{}...".format(int(t)), end=" ")
                    sys.stdout.flush()

            if t + h > b:
                h = b - t

            if False:  # np.any(np.logical_and(x != 0.0, x < self.atol)):
                # do Euler forward step
                dfdt_ = f(x, t)
                m = np.nonzero(dfdt_)
                dt_max = 0.8 * np.min(np.abs(x[m] / dfdt_[m]))

                if self.debug:
                    print("did euler forward step", dt_max)

                x += dt_max * dfdt_
                t += dt_max

                if np.any(np.isnan(x)):
                    raise Exception("New Euler-forward step has nans :'(")
                T = np.append(T, t)
                X = np.append(X, [x], 0)
                continue
                # return ( X, T )

            a2 = 0.5
            b21 = 0.5
            a3 = 0.5
            b31 = 0.0
            b32 = 0.5
            a4 = 1.0
            b41 = 0.0
            b42 = 0.0
            b43 = 1.0
            a5 = 1.0
            b51 = 1.0 / 6.0
            b52 = 1.0 / 3.0
            b53 = 1.0 / 3.0
            b54 = 1.0 / 6.0
            c1_1 = 0.0
            c2_1 = 0.0
            c3_1 = 0.0
            c4_1 = 0.0
            c5_1 = 1.0
            c1_2 = 1.0 / 6.0
            c2_2 = 1.0 / 3.0
            c3_2 = 1.0 / 3.0
            c4_2 = 0.0
            c5_2 = 1.0 / 6.0

            k1 = h * f(x, t)
            k2 = h * f(x + b21 * k1, t + a2 * h)
            k3 = h * f(x + b31 * k1 + b32 * k2, t + a3 * h)
            k4 = h * f(x + b41 * k1 + b42 * k2 + b43 * k3, t + a4 * h)
            k5 = h * f(x + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4, t + a5 * h)

            x_n1 = x + c1_1 * k1 + c2_1 * k2 + c3_1 * k3 + c4_1 * k4 + c5_1 * k5
            x_n2 = x + c1_2 * k1 + c2_2 * k2 + c3_2 * k3 + c4_2 * k4 + c5_2 * k5

            abs_err = np.abs(1.0 / 6.0 * (k4 - k5))

            if self.debug:
                print("dfdz")
                Var.print_formatted(k1 / h)

                # print "error rel:"
                # Var.print_formatted(r/tol)
                print("error abs:")
                Var.print_formatted(abs_err)

                print("sol1")
                # sol1 = x + c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5
                Var.print_formatted(x_n1)

                print("sol2")
                # sol2 =  -(r1 * k1 + r3 * k3 + r4 * k4 + r5 * k5 + r6 * k6 - sol1)
                Var.print_formatted(x_n2)

                print("diff")
                Var.print_formatted(abs_err)

            atol = self.atol
            rtol = self.rtol

            max_total_error = atol + rtol * abs(x)
            try:
                s = (
                    0.84
                    * (np.min((max_total_error[abs_err > 0] / abs_err[abs_err > 0])))
                    ** 0.2
                )
            except ValueError:
                s = 0.5
            h = h * s

            h = min(5.0, h)

            if self.debug:
                Var.print_formatted(max_total_error > abs_err)

            if np.any(np.isnan(k5)) or np.any(np.isnan(k4)):
                # if self.debug:
                print("Found nan, scaling down")
                h *= 0.8

                if h < self.min_step:
                    raise Exception("Step became too small")
            else:

                # if len( np.shape( r ) ) > 0:
                # r = max( r )
                # Var.print_formatted(np.abs(r) <= tol)
                if np.all(abs_err <= max_total_error):
                    t = t + h
                    # x = x + c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5
                    x = x_n1
                    T = np.append(T, t)
                    X = np.append(X, [x], 0)
                    S = np.append(S, x)

        return (X, T)
