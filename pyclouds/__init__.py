import numpy as np

# r, w, T, q_v, q_r, q_l, q_i
class Var:
    """
    Representation of a state vector at a specific altitude
    """

    r = 0
    w = 1
    T = 2
    q_v = 3
    q_r = 4
    q_l = 5
    q_i = 6
    z = 7
    p = 8
    q_pr = 9  # TODO: should this really be a m_precip/m_total?

    names = ["r", "w", "T", "q_v", "q_r", "q_l", "q_i", "z", "p", "q_pr"]
    NUM = len(names)

    @staticmethod
    def print_formatted(v, formatting="%g"):
        print(
            ",\t".join(
                [("%s=" + formatting) % (Var.names[i], v[i]) for i in range(Var.NUM)]
            )
        )

    @staticmethod
    def repr(v, formatting="%g", skip=[]):
        units = {
            "w": "m/s",
            "r": "m",
            "T": "K",
            "z": "m",
            "p": "Pa",
        }
        return ", ".join(
            [
                r"$%s=%g%s$" % (Var.names[i], v[i], units.get(Var.names[i], ""))
                for i in range(Var.NUM)
                if not Var.names[i] in skip
            ]
        )

    @staticmethod
    def make_state(**kwargs):
        s = np.zeros((Var.NUM))
        for k, v in list(kwargs.items()):
            s[getattr(Var, k)] = v

        return s


REQUIRED_POSITIVE = np.zeros((Var.NUM))
# REQUIRED_POSITIVE[Var.r] = 1
# REQUIRED_POSITIVE[Var.q_v] = 1
REQUIRED_POSITIVE[Var.q_l] = 1
# REQUIRED_POSITIVE[Var.q_r] = 1
# REQUIRED_POSITIVE[Var.q_i] = 1


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

from . import plot
