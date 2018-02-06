nz = 200
nvars = 10

# initial conditions
T_base, q_base, r_base, w_base = 300., 0.001, 100., 0.1
z_base = 800.

# ambient conditions


z_profile = np.linspace(0., 100., nz)

F_profile = np.array((nz, nvars))

F0 = (T_base, q_base, r_base, w_base)

F_profile[0] = F0
z = z_base
for n in range(1, nz):
    F_profile[n] = integrate(dFdz, y0=F_profile[n-1], x0=z, x_end=z_profile[n])
