import matplotlib.pyplot as plot

import utils
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

def hydrometeor_profile_plot(F, z, Te, p_e):
    p = p_e
    r = F[:,Var.r]
    w = F[:,Var.w]
    T = F[:,Var.T]
    q_v = F[:,Var.q_v]
    q_l = F[:,Var.q_l]
    q_r = F[:,Var.q_r]

    T = [float(t_) for t_ in T]

    q_v__sat = utils.qv_sat(T=T, p=p)

    plot.figure(figsize=(16,6))
    plot.subplot(131)
    plot.plot(q_v, z, marker='x', label='in-cloud')
    plot.plot(q_v__sat, z, marker='x', label='saturation (in-cloud temp)')
    plot.ylabel('height [m]')
    plot.xlabel('water vapor specific concentration [kg/kg]')
    plot.legend()
    plot.grid(True)
    plot.subplot(132)
    plot.plot(q_l, z, marker='x', label='in-cloud')
    plot.ylabel('height [m]')
    plot.xlabel('liquid water specific concentration [kg/kg]')
    plot.grid(True)
    plot.subplot(133)
    plot.plot(q_r, z, marker='x', label='in-cloud')
    plot.ylabel('height [m]')
    plot.xlabel('rain water specific concentration [kg/kg]')
    plot.grid(True)
