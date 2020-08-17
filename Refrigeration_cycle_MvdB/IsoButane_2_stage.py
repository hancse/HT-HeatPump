import CoolProp
import matplotlib.pyplot as plt
import math
from CoolProp.Plots import StateContainer

<<<<<<< HEAD
refrigerant = 'R600a'

=======
>>>>>>> 2594c847dd446fb39b0c7b74ac936497e42dfb7c
# setting the gasses
HEOS = CoolProp.AbstractState('HEOS', 'IsoButane')


# Setting the variables
evaporation_temp = 0 + 273.15
condensation_temp = 75 + 273.15
superheat_evap = 10
superheat_intercooler = 10

subcooling_intercooler = 5
isentropic_eff = 0.7
compressor_loss = 0.1
condensor_capacity = 50000

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
<<<<<<< HEAD

P7 = P6 - DP_cond
=======
P7 = P6
>>>>>>> 2594c847dd446fb39b0c7b74ac936497e42dfb7c
HEOS.specify_phase(CoolProp.iphase_gas)
HEOS.update(CoolProp.PQ_INPUTS, P7, 0)
h7 = HEOS.hmass()
T7 = HEOS.T()
s7 = HEOS.smass

"""----Subcooling point intercooler ----"""
<<<<<<< HEAD

P8a = P7 - DP_sub
=======
P8a = P7
>>>>>>> 2594c847dd446fb39b0c7b74ac936497e42dfb7c
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
<<<<<<< HEAD

P10a = P1a + DP_int
=======
P10a = P1a
>>>>>>> 2594c847dd446fb39b0c7b74ac936497e42dfb7c
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
<<<<<<< HEAD


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


check_massflow = ((M-m_1)/M == m_2/M == X_2)
Compressor_power_lower_stage = (m_2*(h4a-h2))
Compressor_power_higher_stage = (M*(h4-h2a))
Total_compressor_power = Compressor_power_lower_stage + \
                         Compressor_power_higher_stage
COP = condensor_capacity/Total_compressor_power

print("Total massflow:", round(M, 5), "[kg/s]")
print("Massflow evaporator:", round(m_2, 5), "[kg/s]")
print("Massflow intermediate stage:", round(m_1, 5), "[kg/s]")
print("Fraction of massflow through evaporator:", round(X_2*100, 3), "%")
print("Power lower stage:", round(Compressor_power_lower_stage, 0), "[W]")
print("Power upper stage:", round(Compressor_power_higher_stage, 0), "[W]")
print("COP:", round(COP, 4))
=======
check_massflow = ((M-m_1)/M == m_2/M == X_2)
COP = condensor_capacity/((m_2*(h4a-h2)) + (M*(h4-h2a)))
print("Total massflow:", M)
print("Massflow evaporator:", m_2)
print("Massflow intermediate stage:", m_1)
print("Fraction of massflow through evaporator:", X_2)
print("COP:", COP)
>>>>>>> 2594c847dd446fb39b0c7b74ac936497e42dfb7c

HEOS.build_phase_envelope("dummy")
PE = HEOS.get_phase_envelope_data()
PELabel = 'Propane, x = '
plt.plot([x / HEOS.molar_mass() for x in PE.hmolar_vap],PE.p, '-', label=PELabel)


h=[h1,h2,h4a,h5a,h2a, h4, h6, h7, h8, h10, h1,h2,h4a,h2a, h4, h6, h7, h8a, h10a,h2a]
P=[P1,P2,P4a,P5a,P2a, P4, P6, P7, P8, P10, P1,P2,P4a,P2a, P4, P6, P7, P8a, P10a,P2a]
"""
Total_massflow_calculated = condensor_capacity/(h4hot-h4bub)


deltahsub1 = h4bub-h4sub1
qdotsub1 = deltahsub1*Total_massflow_calculated

deltahintermediate = h3dintermediate - h6
massflow_intermediate = qdotsub1/deltahintermediate

massflow_low = Total_massflow_calculated-massflow_intermediate


power_lower_stage = (h3disintermediate-h2)*massflow_low
power_high_stage = (h4dis-h3)*Total_massflow_calculated

Calculated_COP= condensor_capacity/(power_lower_stage+power_high_stage)
compressor_power = power_lower_stage+power_high_stage
print("Massflow:",Total_massflow_calculated)
print("Calculated_COP:",Calculated_COP)
print("Calculated_Compressor power:",compressor_power)
"""
plt.plot(h,P,'red')
plt.xlabel('Enthalpy [J/kg]')
plt.ylabel('Pressure [kPa]')
plt.yscale('log')
plt.ylim(9e4,5e6)
plt.xlim(320000,700000)
plt.title('Refrigerant cycle for R290/R600 Mixtures')
plt.legend(loc='lower right', shadow=True)

"""----------------------------------------------------"""

cycle_states = StateContainer()

cycle_states[" 1 " ,'H'] = h1
cycle_states[" 1 "]['S'] = s1
cycle_states[" 1 "][CoolProp.iP] = P1
cycle_states[" 1 ",CoolProp.iT] = T1

cycle_states[" 2 " ,'H'] = h2
cycle_states[" 2 "]['S'] = s2
cycle_states[" 2 "][CoolProp.iP] = P2
cycle_states[" 2 ",CoolProp.iT] = T2

cycle_states[" 3a " ,'H'] = h3a
cycle_states[" 3a "]['S'] = s3a
cycle_states[" 3a "][CoolProp.iP] = P3a
cycle_states[" 3a ",CoolProp.iT] = T3a

cycle_states[" 4a " ,'H'] = h4a
cycle_states[" 4a "]['S'] = s4a
cycle_states[" 4a "][CoolProp.iP] = P4a
cycle_states[" 4a ",CoolProp.iT] = T4a

cycle_states[" 5a " ,'H'] = h5a
cycle_states[" 5a "]['S'] = s5a
cycle_states[" 5a "][CoolProp.iP] = P5a
cycle_states[" 5a ",CoolProp.iT] = T5a

cycle_states[" 2a " ,'H'] = h2a
cycle_states[" 2a "]['S'] = s2a
cycle_states[" 2a "][CoolProp.iP] = P2a
cycle_states[" 2a ",CoolProp.iT] = T2a

cycle_states[" 3 " ,'H'] = h3
cycle_states[" 3 "]['S'] = s3
cycle_states[" 3 "][CoolProp.iP] = P3
cycle_states[" 3 ",CoolProp.iT] = T3

cycle_states[" 4 " ,'H'] = h4
cycle_states[" 4 "]['S'] = s4
cycle_states[" 4 "][CoolProp.iP] = P4
cycle_states[" 4 ",CoolProp.iT] = T4

cycle_states[" 5 " ,'H'] = h5
cycle_states[" 5 "]['S'] = s5
cycle_states[" 5 "][CoolProp.iP] = P5
cycle_states[" 5 ",CoolProp.iT] = T5

cycle_states[" 6 " ,'H'] = h6
cycle_states[" 6 "]['S'] = s6
cycle_states[" 6 "][CoolProp.iP] = P6
cycle_states[" 6 ",CoolProp.iT] = T6

cycle_states[" 7 " ,'H'] = h7
cycle_states[" 7 "]['S'] = s7
cycle_states[" 7 "][CoolProp.iP] = P7
cycle_states[" 7 ",CoolProp.iT] = T7

cycle_states[" 8a " ,'H'] = h8a
cycle_states[" 8a "]['S'] = s8a
cycle_states[" 8a "][CoolProp.iP] = P8a
cycle_states[" 8a ",CoolProp.iT] = T8a


cycle_states[" 10a " ,'H'] = h10a
cycle_states[" 10a "]['S'] = s10a
cycle_states[" 10a "][CoolProp.iP] = P10a
cycle_states[" 10a ",CoolProp.iT] = T10a

cycle_states[" 8 " ,'H'] = h8
cycle_states[" 8 "]['S'] = s8
cycle_states[" 8 "][CoolProp.iP] = P8
cycle_states[" 8 ",CoolProp.iT] = T8

cycle_states[" 10 " ,'H'] = h10
cycle_states[" 10 "]['S'] = s10
cycle_states[" 10 "][CoolProp.iP] = P10
cycle_states[" 10 ",CoolProp.iT] = T10

cycle_states[" 11a ", 'H'] = h11a
cycle_states[" 11a "]['S'] = s11a
cycle_states[" 11a "][CoolProp.iP] = P11a
cycle_states[" 11a ", CoolProp.iT] = T11a

print(cycle_states)
