from .. import AttrDict
from . import constants as reference_constants

import numpy as np
import warnings


class BaseParameterisation(object):
    def __init__(self, constants=None):
        if constants is None:
            constants = self.default_constants

        self.constants = AttrDict(constants)


class SaturationVapourPressure(BaseParameterisation):
    default_constants = {
        "p0vs": 611.2,  # [Pa]
        "a0_lq": 17.67,
        "a1_lq": -32.19,
        "a0_ice": 22.587,
        "a1_ice": 0.7,
        "R_d": 287.05,
        "R_v": 461.51,
    }

    def pv_sat_liquid(self, T):
        p0vs = self.constants.p0vs
        a0_lq = self.constants.a0_lq
        a1_lq = self.constants.a1_lq

        return p0vs * np.exp((a0_lq * (T - 273.15) / (T + a1_lq)))

    def pv_sat_ice(self, T):
        p0vs = self.constants.p0vs
        a0_ice = self.constants.a0_ice
        a1_ice = self.constants.a1_ice

        return p0vs * np.exp((a0_ice * (T - 273.15) / (T + a1_ice)))

    def pv_sat(self, T):
        v = np.zeros(np.array(T).shape)

        if v.shape == ():
            if T > 273.15:
                return self.pv_sat_liquid(T)
            else:
                return self.pv_sat_ice(T)
        else:
            T = np.array(T)
            idx_liquid = np.array(T) > 273.15
            idx_ice = idx_liquid == False
            v[idx_liquid] = self.pv_sat_liquid(T[idx_liquid])
            v[idx_ice] = self.pv_sat_ice(T[idx_ice])
            return v

    def qv_sat(self, T, p):
        R_v = self.constants.R_v
        R_d = self.constants.R_d

        pv_sat = self.pv_sat(T=T)
        epsilon = R_d / R_v
        qv_sat = (epsilon * pv_sat) / (p - (1.0 - epsilon) * pv_sat)

        return qv_sat

    def dpsat_dT(self, T):
        if T < 273.15:
            A, B = self.constants.a0_ice, self.constants.a1_ice
        else:
            A, B = self.constants.a0_lq, self.constants.a1_lq

        return self.pv_sat(T=T) * A * (273.15 - B) / ((T - B) ** 2.0)

    def dqv_sat__dT(self, p, T):
        dpsat_dT = self.dpsat_dT(T=T)

        R_v = self.constants.R_v
        R_d = self.constants.R_d

        pv_sat = self.pv_sat(T=T)

        return R_d / R_v * p * dpsat_dT / ((p - pv_sat) ** 2.0)

    def __call__(self, T):
        return self.pv_sat(T)

    def __str__(self):
        constants_label_s = {
            "ATHAM": reference_constants.ATHAM_constants,
            "pyclouds": reference_constants.default_constants,
            "CCFM": reference_constants.CCFM_constants,
        }

        constants_label = ""
        for k, v in list(constants_label_s.items()):
            if self.constants == v:
                constants_label = k + ": "
                break

        constants_label += ", ".join(
            ["%s=%.03g" % (k, v) for (k, v) in list(self.constants.items())]
        )

        return constants_label


class SaturatedAdiabaticLapseRate(BaseParameterisation):
    """
    Compute saturated moist lapse rate, equation taken directly from wikipedia

    https://en.wikipedia.org/wiki/Lapse_rate#Moist_adiabatic_lapse_rate
    """

    default_constants = reference_constants.default_constants

    def __call__(self, p, T):
        pv_sat_ = pv_sat.pv_sat(T=T)
        g = self.constants.get("g")
        L_v = self.constants.get("L_v")
        R_v = self.constants.get("R_v")
        R_d = self.constants.get("R_d")
        eps = R_d / R_v
        cp_d = self.constants.get("cp_d")

        r = eps * pv_sat_ / (p - pv_sat_)

        dTdz_moist = (
            g
            * (1.0 + L_v * r / (R_d * T))
            / (cp_d + L_v ** 2.0 * r * eps / (R_d * T ** 2.0))
        )

        return dTdz_moist


class DynamicViscosity(BaseParameterisation):
    default_constants = {}

    class Implementations:
        ROGERS_AND_YAU = 0
        G_THOMPSON = 1

    def __init__(self, implementation=Implementations.ROGERS_AND_YAU, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        assert implementation in list(
            DynamicViscosity.Implementations.__dict__.values()
        )

        self.implementation = implementation

    @np.vectorize
    def __thompson(T):
        """
        from G. Thompson '07 microphysics scheme
        """
        Tc = T - 273.15

        if Tc > 0.0:
            return (1.718 + 0.0049 * Tc) * 1.0e-5
        else:
            return (1.718 + 0.0049 * Tc - 1.2e-5 * Tc ** 2.0) * 1.0e-5

    def __call__(self, T):
        """
        Sources:
            Rogers & Yau 1989
        """
        if self.implementation == self.Implementations.ROGERS_AND_YAU:
            return 1.72e-5 * (393.0 / (T + 120.0)) * (T / 273.0) ** (3.0 / 2.0)
        elif self.implementation == self.Implementations.G_THOMPSON:
            return self.__thompson(T)
        else:
            raise NotImplementedError

    def __str__(self):
        if self.implementation == self.Implementations.ROGERS_AND_YAU:
            return "Rogers & Yau 1989"
        elif self.implementation == self.Implementations.G_THOMPSON:
            return "G. Thompson '07 microphysics"
        else:
            raise NotImplementedError


class ThermalConductivityCoefficient(BaseParameterisation):
    default_constants = {
        "a_K": 8.0e-5,
        "b_K": 2.4e-2,
    }

    def __call__(self, T):
        return self.constants.a_K * (T - 273.15) + self.constants.b_K

    def __str__(self):

        return "Linear model (%s)" % ", ".join(
            ["%s=%.3g" % (k, v) for (k, v) in list(self.constants.items())]
        )


class WaterVapourDiffusionCoefficient(BaseParameterisation):
    ATHAM_constants = {
        "a": 2.11e-5,
        "b": 1.94,
    }

    # these constants were obtained by fitting against data in Rogers & Yau
    default_constants = {
        "a": 2.20e-5,
        "b": 1.92,
    }

    def __call__(self, T, p):
        a = self.constants.a
        b = self.constants.b
        T0 = 273.15

        # tabulated values in Rogers & Yau are given for p=100kPa=100000Pa reference pressure,
        # have to scale by pressure get the correct diffusivity
        p0 = 100e3

        return a * (T / T0) ** b * p0 / p

    def __str__(self):
        constants_label_s = {
            "ATHAM": self.ATHAM_constants,
            "fitted to Rogers & Yau": self.default_constants,
        }

        constants_label = ""
        for k, v in list(constants_label_s.items()):
            if self.constants == v:
                constants_label = k + ": "
                break

        constants_label += ", ".join(
            ["%s=%.03g" % (k, v) for (k, v) in list(self.constants.items())]
        )

        return "power-law (%s)" % constants_label


class ParametersationsWithSpecificConstants:
    def __wrap(self, parameterisation, constants, subset):
        default_constants = parameterisation.default_constants

        new_constants = {}
        for c_name in list(default_constants.keys()):
            try:
                new_value = constants[subset][c_name]
            except KeyError:
                new_value = constants.get(c_name, None)

            if new_value is None:
                new_value = default_constants.get(c_name)

            new_constants[c_name] = new_value

        return parameterisation(new_constants)

    def __init__(self, constants):
        constants = reference_constants.make_related_constants(constants)
        self.pv_sat = self.__wrap(SaturationVapourPressure, constants, "pv_sat")
        self.Ka = self.__wrap(ThermalConductivityCoefficient, constants, "Ka")
        self.Dv = WaterVapourDiffusionCoefficient()
        self.dTdz_moist = SaturatedAdiabaticLapseRate(constants)
        self.dyn_visc = DynamicViscosity(constants=constants)


pv_sat = SaturationVapourPressure()
Ka = ThermalConductivityCoefficient()
Dv = WaterVapourDiffusionCoefficient()
dTdz = SaturatedAdiabaticLapseRate()
dyn_visc = DynamicViscosity()
