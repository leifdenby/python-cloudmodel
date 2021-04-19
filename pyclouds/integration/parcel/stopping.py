from ... import Var


def very_high_temperature(z, F):
    return F[Var.T] - 320.0


very_high_temperature.direction = 1.0


def parcel_stopped_rising(z, F):
    return F[Var.w]


parcel_stopped_rising.direction = -1.0


def very_large_radius(z, F):
    return F[Var.r] - 10e3


very_large_radius.direction = 1.0


def radius_below_zero(z, F):
    return F[Var.r]


radius_below_zero.direction = -1.0


def height_unphysical(z, F):
    return z - 40e3


height_unphysical.direction = 1.0


def height_below_ground(z, F):
    return z - 1.0


height_below_ground.direction = -1.0


DEFAULT_STOPPING_FUNCTIONS = [
    very_high_temperature,
    parcel_stopped_rising,
    very_large_radius,
    radius_below_zero,
    height_unphysical,
    height_below_ground,
]
