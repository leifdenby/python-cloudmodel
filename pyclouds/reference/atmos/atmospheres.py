from math import sqrt

import reference.atmospheric_flow.gas_properties


class StandardEarthAtmosphere:
    rho = 1.205
    p = 101325.0

    def __init__(self):
        self.gas_properties = reference.atmospheric_flow.gas_properties.AtmosphericAir()

    def getSoundSpeed(self):
        return sqrt(self.p * self.gas_properties.gamma() / self.rho)

    def __str__(self):
        return "$\\rho$=%g $kg/m^3$, p=%g" % (self.rho, self.p)
