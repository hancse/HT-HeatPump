import CoolProp
import matplotlib.pyplot as plt
import math 

#setting the gasses
propane_mass_fraction = 0.5
HEOS = CoolProp.AbstractState('HEOS', 'IsoButane')
#HEOS.set_mass_fractions([0.3, 0.7])

#Setting the points
evaporation_temp = 0 + 273.15                
condensation_temp = 75 + 273.15
superheat_evap = 5
superheat_intercooler = 2 

subcooling_intercooler = 10                              
isentropic_eff = 0.7
compressor_loss = 0.1
condensor_capacity = 14770


"""----Evaporation----"""

T1 = evaporation_temp
HEOS.update(CoolProp.QT_INPUTS, 1, T1)
P1 = HEOS.p()
h1 = HEOS.hmass()
s1 = HEOS.smass()

"""----Superheating----"""

P2 = P1
T2 = superheat_evap + evaporation_temp
HEOS.update(CoolProp.PT_INPUTS, P2, T2)
h2= HEOS.hmass()
s2= HEOS.smass()


"""----100% isentopic compression lower stage----"""
T6 = condensation_temp
HEOS.update(CoolProp.QT_INPUTS, 1, T6)
P6 = HEOS.p()
P1a = math.sqrt(P1*P6)

#create inital guesses for T3s and s3s
T3a = condensation_temp
HEOS.update(CoolProp.PT_INPUTS, P1a, T3a)
s3a = HEOS.smass()

#Using a loop to find real values for H3s and T3s
while (not(s3a < (s2+1)  and s3a > (s2-1))):
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

while (not(h2dummy < (h4a+1000)  and h2dummy > (h4a-1000))):
    if h2dummy < h4a:
        T4a = T4a + 0.1
    else:
        T4a = T4a -0.1
    HEOS.update(CoolProp.PT_INPUTS, P4a, T4a)
    h2dummy = HEOS.hmass()
T4a = HEOS.T()
s4a = HEOS.smass()

"""----Hot point lower stage----"""
h5a = h4a -((h4a-h2)*compressor_loss)
P5a = P1a
h2dummy = 1000
T5a = T4a
while (not(h2dummy < (h5a+1000)  and h2dummy > (h5a-1000))):
    if h2dummy < h5a:
        T5a = T5a + 0.1
    else:
        T5a = T5a - 0.1
    HEOS.update(CoolProp.PT_INPUTS, P2, T5a)
    h2dummy = HEOS.hmass()
T5a = HEOS.T()
s5a = HEOS.smass()

"""---Intermediate boiling point----"""
HEOS.update(CoolProp.PQ_INPUTS, P1a, 1)
T1a = HEOS.T()
#HEOS.update(CoolProp.PT_INPUTS, P1a, T1a)
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
T3 = condensation_temp -10
P3 = P6

#create inital guesses for T3s and s3s
HEOS.update(CoolProp.PT_INPUTS, P3, T3)
s3 = HEOS.smass()

#Using a loop to find real values for H3s and T3s
while (not(s3 < (s2a+1)  and s3 > (s2a-1))):
    if s3 < s2a:
        T3 = T3 + 0.1
    else:
        T3 = T3 -0.1
    HEOS.update(CoolProp.PT_INPUTS, P3, T3)
    s3 = HEOS.smass()
h3 = HEOS.hmass()

"""----Discharge point upper stage----"""
h4 = ((h3-h2a)/isentropic_eff) + h2a
P4 = P6
HEOS.update(CoolProp.PT_INPUTS, P4, T6)
T4 = HEOS.T()
h4test = HEOS.hmass()

while (not(h4test < (h4+1000)  and h4test > (h4-1000))):
    if h4test < h4:
        T4 = T4 + 0.1
    else:
        T4 = T4 -0.1
    HEOS.update(CoolProp.PT_INPUTS, P4, T4)
    h4test = HEOS.hmass()
T4 = HEOS.T()
s4 = HEOS.smass()

"""----Hot point upper stage----"""
h5 = h4-((h4-h2a)*compressor_loss)
P5 = P4
h5test = 1000
T5 = T4
while (not(h5test < (h5+1000)  and h5test > (h5-1000))):
    if h5test < h5:
        T5 = T5 + 0.1
    else:
        T5 = T5 - 0.1
    HEOS.update(CoolProp.PT_INPUTS, P5, T5)
    h5test = HEOS.hmass()
T5 = HEOS.T()
s4hot = HEOS.smass()

"""---Condensation Point---"""
HEOS.update(CoolProp.PQ_INPUTS, P6, 1)
h6 = HEOS.hmass()
s6 = HEOS.smass()

"""---Bubble point---"""
P7 = P6
HEOS.update(CoolProp.PQ_INPUTS, P7, 0)
h7 = HEOS.hmass()
T7 = HEOS.T()
s7 = HEOS.smass

"""----Subcooling point IHE Intermediate----"""
P8a = P7
HEOS.specify_phase(CoolProp.iphase_liquid)
HEOS.update(CoolProp.PQ_INPUTS, P8a, 0)
T8a = T7 - subcooling_intercooler
HEOS.update(CoolProp.PT_INPUTS, P3, T8a)
h8a = HEOS.hmass()
s8a = HEOS.smass()

"""----Subcooling point IHE Intermediate----"""
P8 = P8a
h8 = h8a - (h2-h1)
h8test = 1000
T8 = T8a
while (not(h8test < (h8+1000)  and h8test > (h8-1000))):
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
P10 = P1
HEOS.update(CoolProp.PQ_INPUTS, P10, 1)
HEOS.specify_phase(CoolProp.iphase_gas)
hg = HEOS.hmass()
HEOS.update(CoolProp.PQ_INPUTS, P10, 0)
HEOS.specify_phase(CoolProp.iphase_liquid)
hf = HEOS.hmass()
X10=(h10-hf)/(hg-hf) 

HEOS.specify_phase(CoolProp.iphase_twophase)
HEOS.update(CoolProp.PQ_INPUTS, P10, X10)
T10 = HEOS.T()
s10 = HEOS.smass()

"""---Intermediate expansion part---"""
P10a = P1a
h10a = h8a
HEOS.update(CoolProp.PQ_INPUTS, P10a, 1)
HEOS.specify_phase(CoolProp.iphase_gas)
hg10a = HEOS.hmass()
HEOS.update(CoolProp.PQ_INPUTS, P10a, 0)
HEOS.specify_phase(CoolProp.iphase_liquid)
hf10a = HEOS.hmass()
X10a=(h10a-hf10a)/(hg10a-hf10a)

HEOS.specify_phase(CoolProp.iphase_twophase)
HEOS.update(CoolProp.PQ_INPUTS, P10a, X10a)
T10a = HEOS.T()
s10a = HEOS.smass()

HEOS.build_phase_envelope("dummy")
PE = HEOS.get_phase_envelope_data()
PELabel = 'Propane, x = ' + str(propane_mass_fraction)
plt.plot([x / HEOS.molar_mass() for x in PE.hmolar_vap],PE.p, '-', label=PELabel)


h=[h1,h2,h4a,h2a, h4, h6, h7, h8, h10, h1,h2,h4a,h2a, h4, h6, h7, h8a, h10a,h2a]
P=[P1,P2,P4a,P2a, P4, P6, P7, P8, P10, P1,P2,P4a,P2a, P4, P6, P7, P8a, P10a,P2a]
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
plt.ylim(5e4,10e6)
plt.xlim(200000,750000)
plt.title('Refrigerant cycle for R290/R600 Mixtures')
plt.legend(loc='lower right', shadow=True)
