
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

