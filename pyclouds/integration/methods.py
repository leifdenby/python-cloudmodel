import sys
import numpy as np

from .common import Var

class NewSolver:
    def __init__(self, dFdz, abs_tol, rel_tol, min_step):
        self.dFdz = dFdz
        self.abs_tol = abs_tol
        self.rel_tol = rel_tol
        self.min_step = min_step
        self.debug = False

    def set_initial_condition(self, F0):
        self.F0 = F0

    def solve(self, z, stopping_func):
        f=self.dFdz
        x0=self.F0
        a=z.min()
        b=z.max()
        hmax=10
        hmin=self.min_step


        t = a
        x = x0
        h = 1.0

        # Initialize arrays that will be returned

        T = np.array( [t] )
        X = np.array( [x] )

        S = np.array([np.nan,])

        kk = 0

        while t < b:
            kk += 1
            # Adjust step size when we get to last interval

            if stopping_func(X, T, len(T)-1):
                return (X, T)

            if self.debug:
                print()
                print("kk, t, h", kk, t, h)

                print("start")
                Var.print_formatted(x)
            else:
                if int(t) % 100 == 0:
                    print("{}...".format(int(t)), end=' ')
                    sys.stdout.flush()

            if t + h > b:
                h = b - t;

            if False: #np.any(np.logical_and(x != 0.0, x < self.abs_tol)):
                # do Euler forward step
                dfdt_ = f(x,t)
                m = np.nonzero(dfdt_)
                dt_max = 0.8*np.min(np.abs(x[m]/dfdt_[m]))

                if self.debug:
                    print("did euler forward step", dt_max)

                x += dt_max*dfdt_
                t += dt_max

                if np.any(np.isnan(x)):
                    raise Exception("New Euler-forward step has nans :'(")
                T = np.append( T, t )
                X = np.append( X, [x], 0 )
                continue
                #return ( X, T )

            a2   = 0.5
            b21  = 0.5
            a3   = 0.5
            b31  = 0.0
            b32  = 0.5
            a4   = 1.0
            b41  = 0.0
            b42  = 0.0
            b43  = 1.0
            a5   = 1.0
            b51  = 1./6.
            b52  = 1./3.
            b53  = 1./3.
            b54  = 1./6.
            c1_1 = 0.0
            c2_1 = 0.0
            c3_1 = 0.0
            c4_1 = 0.0
            c5_1 = 1.
            c1_2 = 1./6.
            c2_2 = 1./3.
            c3_2 = 1./3.
            c4_2 = 0.
            c5_2 = 1./6.


            k1 = h * f( x, t )
            k2 = h * f( x + b21 * k1, t + a2 * h )
            k3 = h * f( x + b31 * k1 + b32 * k2, t + a3 * h )
            k4 = h * f( x + b41 * k1 + b42 * k2 + b43 * k3, t + a4 * h )
            k5 = h * f( x + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4, t + a5 * h )


            x_n1 = x + c1_1*k1 + c2_1*k2 + c3_1*k3 + c4_1*k4 + c5_1*k5
            x_n2 = x + c1_2*k1 + c2_2*k2 + c3_2*k3 + c4_2*k4 + c5_2*k5

            abs_err = np.abs(1./6.*(k4 - k5))

            if self.debug:
                print("dfdz")
                Var.print_formatted(k1/h)

                #print "error rel:"
                #Var.print_formatted(r/tol)
                print("error abs:")
                Var.print_formatted(abs_err)

                print("sol1")
                #sol1 = x + c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5
                Var.print_formatted(x_n1)
                
                print("sol2")
                #sol2 =  -(r1 * k1 + r3 * k3 + r4 * k4 + r5 * k5 + r6 * k6 - sol1)
                Var.print_formatted(x_n2)

                print("diff")
                Var.print_formatted(abs_err)

            abs_tol = self.abs_tol
            rel_tol = self.rel_tol

            max_total_error = (abs_tol + rel_tol*abs(x))
            try:
                s = 0.84*(np.min((max_total_error[abs_err > 0]/abs_err[abs_err > 0])))**0.2
            except ValueError:
                s = 0.5
            h = h * s

            h = min(5.0, h)

            if self.debug:
                Var.print_formatted(max_total_error > abs_err)

            if np.any(np.isnan(k5)) or np.any(np.isnan(k4)):
                #if self.debug:
                print("Found nan, scaling down")
                h *= 0.8

                if h < self.min_step:
                    raise Exception("Step became too small")
            else:

                #if len( np.shape( r ) ) > 0:
                    #r = max( r )
                #Var.print_formatted(np.abs(r) <= tol)
                if np.all(abs_err <= max_total_error):
                    t = t + h
                    #x = x + c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5
                    x = x_n1
                    T = np.append( T, t )
                    X = np.append( X, [x], 0 )
                    S = np.append( S, x)

        return ( X, T )
