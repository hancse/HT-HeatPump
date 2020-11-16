import CoolProp
import math

def refrigeration_cycle(evap_temp, condens_temp, condens_capacity):

    # setting the gasses
    HEOS = CoolProp.AbstractState('HEOS', 'IsoButane')

    # Setting the variables
    evaporation_temp = evap_temp + 273.15
    condensation_temp = condens_temp + 273.15
    superheat_evap = 5
    superheat_intercooler = 5

    subcooling_intercooler = 5
    isentropic_eff = 0.7
    compressor_loss = 0.1
    condensor_capacity = condens_capacity

    DP_sup = 5000
    DP_cond = 5000
    DP_sub = 5000
    DP_evap = 5000
    DP_int = 5000

    """----Evaporation----"""
    T1 = evaporation_temp
    HEOS.update(CoolProp.QT_INPUTS, 1, T1)
    HEOS.specify_phase(CoolProp.iphase_gas)
    P1 = HEOS.p()
    h1 = HEOS.hmass()
    s1 = HEOS.smass()

    """----Superheating----"""

    P2 = P1 - DP_sup
    T2 = superheat_evap + evaporation_temp
    HEOS.update(CoolProp.PT_INPUTS, P2, T2)
    h2 = HEOS.hmass()
    s2 = HEOS.smass()

    """----Intermediate stage pressuer calculation----"""
    T6 = condensation_temp
    HEOS.update(CoolProp.QT_INPUTS, 1, T6)
    P6 = HEOS.p()
    P1a = math.sqrt(P2*P6)

    # Create inital guesses for T3s and s3s
    P3a = P1a
    T3a = condensation_temp
    HEOS.update(CoolProp.PT_INPUTS, P1a, T3a)
    s3a = HEOS.smass()

    # Using a loop to find real values for H3s and T3s
    while (not(s3a < (s2+1) and s3a > (s2-1))):
        if s3a < s2:
            T3a = T3a + 0.1
        else:
            T3a = T3a - 0.1
        HEOS.update(CoolProp.PT_INPUTS, P1a, T3a)
        s3a = HEOS.smass()
    h3a = HEOS.hmass()

    """----Discharge point lower stage----"""
    h4a = ((h3a-h2)/isentropic_eff) + h2
    P4a = P1a
    HEOS.update(CoolProp.PT_INPUTS, P4a, T3a)
    T4a = HEOS.T()
    h2dummy = HEOS.hmass()

    while (not(h2dummy < (h4a + 1000) and h2dummy > (h4a - 1000))):
        if h2dummy < h4a:
            T4a = T4a + 0.1
        else:
            T4a = T4a - 0.1
        HEOS.update(CoolProp.PT_INPUTS, P4a, T4a)
        h2dummy = HEOS.hmass()
    T4a = HEOS.T()
    s4a = HEOS.smass()

    """----Hot point lower stage----"""
    h5a = h4a - ((h4a - h2) * compressor_loss)
    P5a = P1a
    h2dummy = 1000
    T5a = T4a
    while (not(h2dummy < (h5a+1000) and h2dummy > (h5a-1000))):
        if h2dummy < h5a:
            T5a = T5a + 0.1
        else:
            T5a = T5a - 0.1
        HEOS.update(CoolProp.PT_INPUTS, P5a, T5a)
        h2dummy = HEOS.hmass()
    T5a = HEOS.T()
    s5a = HEOS.smass()

    """---Intermediate boiling point----"""
    HEOS.update(CoolProp.PQ_INPUTS, P1a, 1)
    T1a = HEOS.T()
    h1a = HEOS.hmass()
    s1a = HEOS.smass()

    """---Intermediate superheat point----"""
    P2a = P1a
    HEOS.update(CoolProp.PQ_INPUTS, P2a, 1)
    T2a = HEOS.T() + superheat_intercooler
    HEOS.update(CoolProp.PT_INPUTS, P2a, T2a)
    h2a = HEOS.hmass()
    s2a = HEOS.smass()

    """----100% isentopic compression upper stage----"""
    T3 = condensation_temp
    P3 = P6

    # Create inital guesses for T3s and s3s
    HEOS.update(CoolProp.PT_INPUTS, P3, T3)
    s3 = HEOS.smass()

    # Using a loop to find real values for H3s and T3s
    while (not(s3 < (s2a+1) and s3 > (s2a-1))):
        if s3 < s2a:
            T3 = T3 + 0.1
        else:
            T3 = T3 - 0.1
        HEOS.update(CoolProp.PT_INPUTS, P3, T3)
        s3 = HEOS.smass()
    h3 = HEOS.hmass()

    """----Discharge point upper stage----"""
    h4 = ((h3-h2a)/isentropic_eff) + h2a
    P4 = P6
    HEOS.update(CoolProp.PT_INPUTS, P4, T6)
    T4 = HEOS.T()
    h4test = HEOS.hmass()

    while (not(h4test < (h4 + 1000) and h4test > (h4 - 1000))):
        if h4test < h4:
            T4 = T4 + 0.1
        else:
            T4 = T4 - 0.1
        HEOS.update(CoolProp.PT_INPUTS, P4, T4)
        h4test = HEOS.hmass()
    T4 = HEOS.T()
    s4 = HEOS.smass()

    """----Hot point upper stage----"""
    h5 = h4-((h4-h2a)*compressor_loss)
    P5 = P4
    h5test = 1000
    T5 = T4
    while (not(h5test < (h5 + 1000) and h5test > (h5 - 1000))):
        if h5test < h5:
            T5 = T5 + 0.1
        else:
            T5 = T5 - 0.1
        HEOS.update(CoolProp.PT_INPUTS, P5, T5)
        h5test = HEOS.hmass()
    T5 = HEOS.T()
    s5 = HEOS.smass()

    """---Condensation Point---"""
    HEOS.update(CoolProp.PQ_INPUTS, P6, 1)
    h6 = HEOS.hmass()
    s6 = HEOS.smass()

    """---Bubble point---"""
    P7 = P6 - DP_cond
    HEOS.specify_phase(CoolProp.iphase_gas)
    HEOS.update(CoolProp.PQ_INPUTS, P7, 0)
    h7 = HEOS.hmass()
    T7 = HEOS.T()
    s7 = HEOS.smass

    """----Subcooling point intercooler ----"""
    P8a = P7 - DP_sub
    HEOS.specify_phase(CoolProp.iphase_liquid)
    HEOS.update(CoolProp.PQ_INPUTS, P8a, 0)
    T8a = T7 - subcooling_intercooler
    HEOS.update(CoolProp.PT_INPUTS, P3, T8a)
    h8a = HEOS.hmass()
    s8a = HEOS.smass()

    """----Subcooling point desuperheater----"""
    P8 = P8a
    h8 = h8a - (h2-h1)
    h8test = 1000
    T8 = T8a
    while (not(h8test < (h8+1000) and h8test > (h8-1000))):
        if h8test < h8:
            T8 = T8 + 0.1
        else:
            T8 = T8 - 0.1
        HEOS.specify_phase(CoolProp.iphase_liquid)
        HEOS.update(CoolProp.PT_INPUTS, P8, T8)
        h8test = HEOS.hmass()
    T8 = HEOS.T()
    s8 = HEOS.smass()

    """----Flash point----"""
    h9 = h8

    """---I point----"""
    h10 = h9
    P10 = P1 + DP_evap
    HEOS.update(CoolProp.PQ_INPUTS, P10, 1)
    HEOS.specify_phase(CoolProp.iphase_gas)
    hg = HEOS.hmass()
    HEOS.update(CoolProp.PQ_INPUTS, P10, 0)
    HEOS.specify_phase(CoolProp.iphase_liquid)
    hf = HEOS.hmass()
    X10 = (h10-hf) / (hg-hf)

    HEOS.specify_phase(CoolProp.iphase_twophase)
    HEOS.update(CoolProp.PQ_INPUTS, P10, X10)
    T10 = HEOS.T()
    s10 = HEOS.smass()

    """---Intermediate expansion part---"""
    P10a = P1a + DP_int
    h10a = h8a
    HEOS.update(CoolProp.PQ_INPUTS, P10a, 1)
    HEOS.specify_phase(CoolProp.iphase_gas)
    hg10a = HEOS.hmass()
    HEOS.update(CoolProp.PQ_INPUTS, P10a, 0)
    HEOS.specify_phase(CoolProp.iphase_liquid)
    hf10a = HEOS.hmass()
    X10a = (h10a-hf10a) / (hg10a-hf10a)

    HEOS.specify_phase(CoolProp.iphase_twophase)
    HEOS.update(CoolProp.PQ_INPUTS, P10a, X10a)
    T10a = HEOS.T()
    s10a = HEOS.smass()

    M = condensor_capacity/(h5-h7)
    X_2 = (h2a-h10a-(h7-h8a))/(h5a-h10a)
    m_2 = M*X_2
    m_1 = M - m_2
    P11 = P10a
    h11 = ((M/m_1)*(h7-h8a))+h10a
    HEOS.update(CoolProp.HmassP_INPUTS, h11, P11)
    s11 = HEOS.smass()
    T11 = HEOS.T()

    """---Intermediate Flash Point---"""
    T11a = T1a
    HEOS.specify_phase(CoolProp.iphase_liquid)
    HEOS.update(CoolProp.QT_INPUTS, 0, T11a)
    P11a = HEOS.p()
    h11a = HEOS.hmass()
    s11a = HEOS.smass()

    Compressor_power_lower_stage = (m_2*(h4a-h2))
    Compressor_power_higher_stage = (M*(h4-h2a))
    Total_compressor_power = Compressor_power_lower_stage + \
                             Compressor_power_higher_stage
    COP = condensor_capacity/Total_compressor_power
    return COP

test = refrigeration_cycle(4, 72, 50000)


print(test)