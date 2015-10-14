
default_constants = {
    "R_d": 287.05,
    "R_v": 461.51,
    "L_v": 2.5008e6,
    "L_s": 2.8345e6,
    "cp_d": 1005.46,
    "cp_v": 1859.0,
    "rho_l": 1000.,
    "rho_i": 500.,
    "g": 9.80665,
    "cv_d": 717.60, # J/kg/K
    "cv_v": 1402.5, # J/kg/K
}

CCFM_constants = {
    "R_d": 287.05,  # 'Rd'
    "R_v": 461.51,  # 'Rv'
    "L_v": 2.5008e6,  # 'aLv'
    "L_s": 2.8345e6,  # 'aLs'
    "cp_d": 1005.46,  # 'Cpd'
    "cp_v": 1859.0,  # not specified in CCFM?
    "rho_l": 1000.,  # 'Rho_w'
    "rho_i": 500.,  # 'Rho_i'
    "g": 9.80665,  # 'g'
}

# r, w, T, q_v, q_r, q_l, q_i
class Var:
    r = 0
    w = 1
    T = 2
    q_v = 3
    q_r = 4
    q_l = 5
    q_i = 6

    names = ['r', 'w', 'T', 'q_v', 'q_r', 'q_l', 'q_i']
    NUM = len(names)

    @staticmethod
    def print_formatted(v):
        print ",\t".join(["%s=%f" % (Var.names[i], v[i]) for i in range(Var.NUM)])

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

