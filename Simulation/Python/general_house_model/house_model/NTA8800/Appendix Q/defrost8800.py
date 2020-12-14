# This is described in appendix Q3.8 of the NTA 8800 on page 897
import numpy as np
import matplotlib
# matplotlib.use('TkAgg')
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


## NOG NIET KLAAR!!! PvK


def defrost8800(T_evap, COP_par, COP_A2W35=None ):
    Tevap_ref = 2 # degrees Celsius
    Tcond_ref = 35 # degrees Celsius
    COP_nofrost = COP_par[0] + COP_par[1] * Tevap_ref + COP_par[2] * Tcond_ref
    if COP_A2W35:
        frost_factor_0 = COP_A2W35 /  COP_nofrost
    else:
        frost_factor_0 = 0.75

    # The correction for frosting Frost_factor( T_air )  now becomes:
    if (T_evap <= -7) or (T_evap >= 7):
        frost_factor = 1
    elif T_evap <= 2:
        frost_factor =  ((T_evap + 7) /9) * frost_factor_0
    else:
        frost_factor = ((7- T_evap) / 5) * frost_factor_0


    # The  COP then becomes:
    # COP_HP(T_air, T_water) = COP_HP_no_frost(T_air, T_water) * Frost factor( t_air)
    # The Power becomes:
    # Power_HP(T_air, T_water)  = Power_HP_no_frost(T_air, T_water)  * Frost factor( t_air)
    # Remark: no separate frost factor for the power is specified. This is an approximation!
    # These measurements are valid for full power

    return frost_factor


if __name__ == "__main__":
    par = np.c_[[5.0, 0.10, -0.05]]
    Tair = np.arange(-10, 11, 1.0)
    # ff = np.zeros(np.shape(Tair))
    ff = np.empty(np.shape(Tair))
    list = []
    for T in Tair:
        list.append(defrost8800(T, par))

    fig = plt.figure()
    plt.plot(Tair, list)
    plt.show()