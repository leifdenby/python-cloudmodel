import numpy as np

default_constants = {
    "R_d": 287.05,
    "R_v": 461.51,
    "L_v": 2.5008e6,
    "L_s": 2.8345e6,
    "cp_d": 1005.46,
    "cv_d": 717.60,  # J/kg/K
    "cp_v": 1859.0,
    "cv_v": 1402.5,  # J/kg/K
    "cp_l": 4183.0,
    "cv_l": 4183.0,  # same as cp as liquid is assumed incompressible
    "cp_i": 2103.0,  # taken from ATHAM
    "rho_l": 1000.0,
    "rho_i": 500.0,
    "g": 9.80665,
    "pv_sat": {  # constants for calculating saturation vapour pressure with Teten's formula
        "p0vs": 611.2,
        "a0_lq": 17.67,
        "a1_lq": -32.19,
    },
}

CCFM_constants = {
    "R_d": 287.05,  # 'Rd'
    "R_v": 461.51,  # 'Rv'
    "L_v": 2.5008e6,  # 'aLv'
    "L_s": 2.8345e6,  # 'aLs'
    "cp_d": 1005.46,  # 'Cpd'
    "cp_v": 1859.0,  # not specified in CCFM?
    "rho_l": 1000.0,  # 'Rho_w'
    "rho_i": 500.0,  # 'Rho_i'
    "g": 9.80665,  # 'g'
    "pv_sat": {
        "p0vs": 610.78,
        "a0_lq": 17.269,
        "a1_lq": -35.86,
        "a0_ice": 21.875,
        "a1_ice": -7.66,
    },
}

ECHAM6_3_constants = {
    "R_d": 287.04,  # 'Rd'
    "R_v": 461.51,  # 'Rv'
    "L_v": 2.5008e6,  # 'aLv'
    "L_s": 2.8345e6,  # 'aLs'
    "cp_d": 1004.64,  # 'Cpd'
    "cp_v": 1869.46,  # not specified in CCFM?
    "cp_i": 2106.0,  # Not used in CCFM?
    "rho_l": 1000.0,  # 'Rho_w'
    "rho_i": 917.0,  # 'Rho_i'
    "g": 9.80665,  # 'g'
    "pv_sat": {
        "p0vs": 610.78,
        "a0_lq": 17.269,
        "a1_lq": -35.86,
        "a0_ice": 21.875,
        "a1_ice": -7.66,
    },
}

ATHAM_constants = {
    "cp_d": 1004.64,
    "cv_d": 717.60,
    "cp_v": 1870.0,
    "cv_v": 1402.5,
    "L_v": 2500800.0,
    "cp_l": 4183.0,
    "cp_i": 2103.0,
    "rho_l": 1000.0,
    "rho_i": 917.0,
    "pv_sat": {
        "p0vs": 610.7,
        "a0_lq": 17.25,
        "a1_lq": -36.0,
        "a0_ice": 22.33,
        "a1_ice": -2.0,
    },
    "g": 9.80616,
}

UCLALES_constants = {
    "cp_d": 1005.0,
    "R_d": 187.04,
    "L_v": 2500000.0,
    "cp_l": 4183.0,
    "cp_i": 2103.0,
    "rho_l": 1000.0,
    "rho_i": 900.0,
    "pv_sat": {
        "p0vs": 610.78000,
        "a0_lq": 17.2693882,
        "a1_lq": -35.860000,
        "a0_ice": 21.8745574,
        "a1_ice": -7.66,
    },
    "g": 9.8,
}


def make_related_constants(constants):
    if "cp_d" in constants and "cp_v" in constants and not "R_d" in constants:
        constants["R_d"] = constants["cp_d"] - constants["cv_d"]
        constants["R_v"] = constants["cp_v"] - constants["cv_v"]
    if "cp_i" in constants and not "cv_i" in constants:
        constants["cv_i"] = constants["cp_i"]
    if "cp_l" in constants and not "cv_l" in constants:
        constants["cv_l"] = constants["cp_l"]
    if all([v in constants for v in ("L_s", "L_v")]) and not "L_f" in constants:
        constants["L_f"] = constants["L_s"] - constants["L_v"]
    return constants


um_constants = {
    "cp_d": 1004.64,
    "cv_d": 717.60,
    "cp_v": 1864.0,
    "cv_v": 1402.55,
    "cp_l": 4183.0,
    "cv_l": 4183.0,
    "cp_i": 2103.0,
    "cv_i": 2103.0,
    "L_v": 2500.8,
    "rho_l": 1000.0,
    "pv_sat": {  # constants for calculating saturation vapour pressure with Teten's formula
        "p0vs": 610.7,
        "a0_lq": 17.25,
        "a1_lq": 36.0,
    },
}

# class MySolver(odespy.RKFehlberg):

# def advance(self):
# f, n, rtol, atol = self.f, self.n, self.rtol, self.atol
# u_n, t_n, t_np1 = self.u[n], self.t[n], self.t[n+1]
# dt = t_np1 - t_n

# dt_max__all = np.abs(f(u_n, t_n)/u_n)
# dt_max__all[REQUIRED_POSITIVE == 0] = 1.0e31

# self.max_step = min(np.min(dt_max__all), self.first_step)
# print(self.max_step, dt_max__all)
# self.first_step = self.max_step
# import ipdb
# with ipdb.launch_ipdb_on_exception():
# return odespy.RKFehlberg.advance(self)
