"""
house model base on 2R2C model with a buffervessel and a radiator
"""

from scipy.integrate import odeint       # ODE solver
import numpy as np                       # linear algebra

def model_buffervessel(x, t, T_outdoor, Q_internal, Q_solar, SP_T, Qinst, CF, Rair_outdoor, Rair_wall, Cair, Cwall, mdot, UAradiator, Crad, Cbuffervessel, cpwater):
    """model function for scipy.integrate.odeint.

    :param x:            (array):   variable array dependent on time with the vairable Air temperature, Wall temperature Return water temperature and buffervessel temperature
    :param t:            (float):
    :param T_outdoor:    (float):  Outdoor temperature in degree C
    :param Q_internal:   (float):  Internal heat gain in [W]
    :param Q_solar:      (float):  Solar irradiation on window [W]
    :param SP_T:         (float):  Setpoint tempearature from thermostat. [C]
    :param Qinst:        (float):  Heating power delivered to the buffervessel [W]
    :param CF:           (float):  factor of Q_solar heat transferred to the air (Unitless)
    :param Rair_outdoor: (float):  Thermal resistance from indoor air to outdoor air [K/W]
    :param Rair_wall:    (float):  Thermal resistance from indoor air to the wall [K/W]
    :param Cair:         (float):  Thermal capacity of the air
    :param Cwall:        (float):  Thermal capacity of the wall
    :param mdot:         (float):  waterflow in the radiator [kg/s]
    :param UAradiator    (float):  Heat transfer coeffiecient of the radiator 
    :return:             (array):  Difference over of the variables in x      

    x,t: ode input function func : callable(x, t, ...) or callable(t, x, ...)
    Computes the derivative of y at t.
    If the signature is ``callable(t, y, ...)``, then the argument tfirst` must be set ``True``.
    """

    # States :

    Tair = x[0]
    Twall = x[1]
    Treturn = x[2]
    Tbuffervessel = x[3]
 
    # Equations :
        
    Tairdt = ((T_outdoor - Tair) / Rair_outdoor + (Twall - Tair) / Rair_wall + UAradiator*(Treturn-Tair) + Q_internal + CF * Q_solar) / Cair
    Twalldt = ((Tair - Twall) / Rair_wall + (1 - CF) * Q_solar) / Cwall
    Treturndt = ((mdot*cpwater*(Tbuffervessel-Treturn)) + UAradiator*(Tair-Treturn)) / Crad
    Tbuffervesseldt = (Qinst + (cpwater*mdot*(Treturn-Tbuffervessel)))/Cbuffervessel

    return [Tairdt, Twalldt, Treturndt, Tbuffervesseldt]


def house_buffervessel(T_outdoor, Q_internal, Q_solar, SP_T, time_sim, CF,
          Rair_outdoor, Rair_wall, Cair, Cwall, mdot, UAradiator, Crad, Cbuffervessel, cpwater):
    """Compute air and wall tempearature inside the house.

    :param T_outdoor:    (array):  Outdoor temperature in degree C
    :param Q_internal:   (array):  Internal heat gain in w.
    :param Q_solar:      (array):  Solar irradiation on window [W]
    :param SP_T:         (array):  Setpoint tempearature from thermostat.
    :param time_sim:     (array)  :  simulation time

    :param CF:
    :param Rair_outdoor:
    :param Rair_wall:
    :param Cair:
    :param Cwall:
    :return:             tuple :  Tuple containing (Tair, Twall):

                - Tair (float):   air temperature inside the house in degree C.
                - Twall (float):  wall temperature inside the house in degree C

    Qinst ?	  (array):  instant heat from heat source such as HP or boiler [W].

    """

    Tair0 = 20
    Twall0 = 20
    Treturn0 = 40
    Tbuffervessel0 = 60

    y0 = [Tair0, Twall0, Treturn0, Tbuffervessel0]

    t = time_sim           # Define Simulation time with sampling time

    Tair = np.ones(len(t)) * Tair0
    Twall = np.ones(len(t)) * Twall0
    Treturn = np.ones(len(t)) * Treturn0
    Tbuffervessel = np.ones(len(t)) * Tbuffervessel0
    consumption = np.ones(len(t))
    kp = 300
    
    for i in range(len(t)-1):

        err = 80 - Tbuffervessel[i]
        err2 = SP_T[i+1] - Tair[i]
        Qinst = err * kp
        Qinst = np.clip(Qinst, 0, 12000)
        
        if (err2 > 0):
            mdot = 0.1
        else:
            mdot = 0

        inputs = (T_outdoor[i], Q_internal[i], Q_solar[i], SP_T[i], Qinst, CF,
                  Rair_outdoor, Rair_wall, Cair, Cwall, mdot, UAradiator, Crad, Cbuffervessel, cpwater)
        ts = [t[i], t[i+1]]
        y = odeint(model_buffervessel, y0, ts, args=inputs)

        Tair[i+1] = y[-1][0]
        Twall[i+1] = y[-1][1]
        Treturn[i+1] = y[-1][2]
        Tbuffervessel[i+1] = y[-1][3]

        # Adjust initial condition for next loop
        y0 = y[-1]
        consumption[i] = Qinst
    return Tair, Twall, Treturn, Tbuffervessel, consumption

