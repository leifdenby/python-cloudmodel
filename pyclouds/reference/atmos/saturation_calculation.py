import numpy as np


epsilon = 0.622

# Teten's formula for saturation vapour pressure (constants from Bechtold's notes)
p0vs = 611.2  # [Pa]
a0_lq = 17.67
a1_lq = -32.19
a0_ice = 22.587
a1_ice = 0.7

# over liquid water
pv_sat_lq = lambda T: p0vs * np.exp((a0_lq * (T - 273.15) / (T + a1_lq)))
# over ice
pv_sat_ice = lambda T: p0vs * np.exp((a0_ice * (T - 273.15) / (T + a1_ice)))


def pv_sat(T):
    v = np.zeros(np.array(T).shape)

    if v.shape == ():
        if T > 273.15:
            return pv_sat_lq(T)
        else:
            return pv_sat_ice(T)
    else:
        T = np.array(T)
        idx_liquid = np.array(T) > 273.15
        idx_ice = idx_liquid == False
        v[idx_liquid] = pv_sat_lq(T[idx_liquid])
        v[idx_ice] = pv_sat_ice(T[idx_ice])
        return v


qv = lambda T, p, pv: (epsilon * pv) / (p - (1 - epsilon) * pv)


def qv_sat(T, p):
    pv = pv_sat(T)
    return qv(T=T, p=p, pv=pv)
