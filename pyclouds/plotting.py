import matplotlib.pyplot as plot

from cloud_equations import Var

def profile_plot(F, z, Te):
    r = F[:,Var.r]
    w = F[:,Var.w]
    T = F[:,Var.T]

    plot.figure(figsize=(16,6))
    plot.subplot(131)
    plot.plot(r, z, marker='x')
    plot.ylabel('height [m]')
    plot.xlabel('radius [m]')
    plot.grid(True)
    plot.subplot(132)
    plot.plot(T, z, marker='x', label='in-cloud')
    plot.plot(Te(z), z, marker='x', label='environment')
    plot.legend()
    plot.xlabel('temperature [K]')
    plot.ylabel('height [m]')
    plot.grid(True)
    plot.subplot(133)
    plot.plot(w, z, marker='x')
    plot.ylabel('height [m]')
    plot.xlabel('vertical velocity [m/s]')
    plot.grid(True)
