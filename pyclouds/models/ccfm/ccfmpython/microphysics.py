from math import exp


Rd = 287.05
Rv = 461.51
cpd = 1005.46
aLv = 2.5008e6
aLs = 2.8345e6
aLf = aLs - aLv
Tnull = 273.15


BASEBUO = 0.5

C1ES = 610.78
C2ES = C1ES * Rd / Rv
C3LES = 17.269
C4LES = 35.86
C3IES = 21.875
C4IES = 7.66

C5LES = C3LES * (Tnull - C4LES)
C5IES = C3IES * (Tnull - C4IES)


def moist_adjust(tem, q_v, prs):
    q_sat = lua(max(180, min(tem, 370))) / prs

    q_sat = min(0.5, q_sat)
    cor = 1 / (1 - (Rv / Rd - 1) * q_sat)
    q_sat = q_sat * cor

    cond = (q_v - q_sat) / (1 + q_sat * cor * lub(max(180, min(tem, 370))))

    cond = max(cond, 0)

    # increase temperature through latent heat of condensation
    tem = tem + luc(tem) * cond

    # remove condensed water from water vapor
    q_v = q_v - cond

    if cond > 0:

        # calculate specific humidity from saturation water vapor pressure using Tetens formula
        q_sat = lua(max(180, min(tem, 370))) / prs

        q_sat = min(0.5, q_sat)
        cor = 1 / (1 - (Rv / Rd - 1) * q_sat)
        q_sat = q_sat * cor

        cond = (q_v - q_sat) / (1 + q_sat * cor * lub(max(180, min(tem, 370))))

        cond = max(cond, 0)

        tem = tem + luc(tem) * cond

        q_v = q_v - cond

    return tem, q_v


def lua(T):
    if T > Tnull:
        CVM3 = C3LES
        CVM4 = C4LES
    else:
        CVM3 = C3IES
        CVM4 = C4IES

    return C2ES * exp(CVM3 * (T - Tnull) / (T - CVM4))


def lub(T):

    if T > Tnull:
        CVM4 = C4LES
        CVM5 = C5LES * aLv / cpd
    else:
        CVM4 = C4IES
        CVM5 = C5IES * aLs / cpd

    return CVM5 / (T - CVM4) ** 2


def luc(T):
    if T > Tnull:
        # latent heat of evaporation / dry air spec. heat cap. at contant pressure
        # J/kg / ( J/kg/K ) = [K*kg/kg]
        luc = aLv / cpd
    else:
        luc = aLs / cpd

    return luc
