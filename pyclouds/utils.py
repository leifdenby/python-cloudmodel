import numpy as np

# Heat capacities
cp_d = 1004.64 # J/kg/K
cv_d = 717.60 # J/kg/K
cp_v = 1874.0 # J/kg/K
cv_v = 1402.5 # J/kg/K

# Teten's formula for saturation vapour pressure (constants from Bechtold's notes)
# over liquid water
p0vs = 611.2  # [Pa]
a0_lq = 17.67
a1_lq = -32.19
a0_ice = 22.587
a1_ice = 0.7

pv_sat_lq = lambda T: p0vs*np.exp((a0_lq*(T-273.15)/(T+a1_lq)))
# over ice
pv_sat_ice = lambda T: p0vs*np.exp((a0_ice*(T-273.15)/(T+a1_ice)))

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

# pv_sat = lambda T: 611.2*np.exp((17.67*(T-273.15)/(T-32.19)))

Rd = cp_d - cv_d
Rv = cp_v - cv_v

epsilon = Rd/Rv

cp = lambda qv: (1-qv)*cp_d + qv*cp_v
cv = lambda qv: (1-qv)*cv_d + qv*cv_v
R = lambda qv: cp(qv) - cv(qv)
kappa = lambda qv: R(qv)/cp(qv)

qv = lambda T, p, pv: (epsilon*pv)/(p-(1-epsilon)*pv)
def qv_sat(T, p):
    pv = pv_sat(T)
    return qv(T=T, p=p, pv=pv)
