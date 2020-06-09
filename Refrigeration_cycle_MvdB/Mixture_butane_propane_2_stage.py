import CoolProp
import matplotlib.pyplot as plt
import math 

#setting the gasses
propane_mass_fraction = 0.5
HEOS = CoolProp.AbstractState('HEOS', 'IsoButane')
#HEOS.set_mass_fractions([0.3, 0.7])

#Setting the points
evaporation_temp = 0 + 273.15                
condensation_temp = 80 + 273.15
superheat = 5                               
isentropic_eff = 0.7
compressor_loss = 0.1
condensor_capacity = 14770


"""----Evaporation----"""

T1 = evaporation_temp
HEOS.update(CoolProp.QT_INPUTS, 1, T1)
HEOS.specify_phase(CoolProp.iphase_gas)
P1 = HEOS.p()
h1 = HEOS.hmass()
s1 = HEOS.smass()

"""----Superheating----"""

P2 = P1
T2 = superheat + evaporation_temp
HEOS.update(CoolProp.PT_INPUTS, P2, T2)
h2= HEOS.hmass()
s2= HEOS.smass()


"""----100% isentopic compression lower stage----"""
T4 = condensation_temp
HEOS.update(CoolProp.QT_INPUTS, 1, T4)
P4 = HEOS.p()
P3sintermediate = math.sqrt(P1*P4)

#create inital guesses for T3s and s3s
T3sintermediate = condensation_temp
HEOS.update(CoolProp.PT_INPUTS, P3sintermediate, T3sintermediate)
s3sintermediate = HEOS.smass()

#Using a loop to find real values for H3s and T3s
while (not(s3sintermediate < (s2+1)  and s3sintermediate > (s2-1))):
    if s3sintermediate < s2:
        T3sintermediate = T3sintermediate + 0.1
    else:
        T3sintermediate = T3sintermediate -0.1
    HEOS.update(CoolProp.PT_INPUTS, P3sintermediate, T3sintermediate)
    s3sintermediate = HEOS.smass()
h3sintermediate = HEOS.hmass()

"""----Discharge point lower stage----"""
h3disintermediate = ((h3sintermediate-h2)/isentropic_eff) + h2
P3disintermediate = P3sintermediate
HEOS.update(CoolProp.PT_INPUTS, P3disintermediate, T3sintermediate)
T3disintermediate = HEOS.T()
h3test = HEOS.hmass()

while (not(h3test < (h3disintermediate+1000)  and h3test > (h3disintermediate-1000))):
    if h3test < h3disintermediate:
        T3disintermediate = T3disintermediate + 0.1
    else:
        T3disintermediate = T3disintermediate -0.1
    HEOS.update(CoolProp.PT_INPUTS, P3disintermediate, T3disintermediate)
    h3test = HEOS.hmass()
T3disintermediate = HEOS.T()
s3disintermediate = HEOS.smass()

"""----Hot point lower stage----"""
h3hotintermediate = h3disintermediate-((h3disintermediate-h2)*compressor_loss)
P3hotintermediate = P3disintermediate
h3test = 1000
T3hotintermediate = T3disintermediate
while (not(h3test < (h3hotintermediate+1000)  and h3test > (h3hotintermediate-1000))):
    if h3test < h3hotintermediate:
        T3hotintermediate = T3hotintermediate + 0.1
    else:
        T3hotintermediate = T3hotintermediate - 0.1
    HEOS.update(CoolProp.PT_INPUTS, P3hotintermediate, T3hotintermediate)
    h3test = HEOS.hmass()
T3hotintermediate = HEOS.T()
s3hotintermediate = HEOS.smass()

"""---Intermediate condensation point----"""
P3dintermediate = P3hotintermediate
HEOS.update(CoolProp.PQ_INPUTS, P3dintermediate, 1)
h3dintermediate = HEOS.hmass()
T3dintermediate = HEOS.T()

"""---Mixing of the 2 stages using IHE---"""
h3 = ((h3hotintermediate-h3dintermediate)/2)+h3dintermediate
P3 = P3dintermediate

h3test = 1000
T3 = T3hotintermediate

while (not(h3test < (h3+1000)  and h3test > (h3-1000))):
    if h3test < h3:
        T3 = T3 + 0.1
    else:
        T3 = T3 - 0.1
    HEOS.update(CoolProp.PT_INPUTS, P3, T3)
    h3test = HEOS.hmass()

T3 = HEOS.T()
s3 = HEOS.smass()

"""----100% isentopic compression upper stage----"""
T4s = condensation_temp
HEOS.update(CoolProp.QT_INPUTS, 1, T4s)
P4s = HEOS.p()

#create inital guesses for T3s and s3s
HEOS.update(CoolProp.PT_INPUTS, P4s, T4s)
s4s = HEOS.smass()

#Using a loop to find real values for H3s and T3s
while (not(s4s < (s3+1)  and s4s > (s3-1))):
    if s4s < s3:
        T4s = T4s + 0.1
    else:
        T4s = T4s -0.1
    HEOS.update(CoolProp.PT_INPUTS, P4s, T4s)
    s4s = HEOS.smass()
h4s = HEOS.hmass()

"""----Discharge point upper stage----"""
h4dis = ((h4s-h3)/isentropic_eff) + h3
P4dis = P4s
HEOS.update(CoolProp.PT_INPUTS, P4s, T4s)
T4dis = HEOS.T()
h4test = HEOS.hmass()

while (not(h4test < (h4dis+1000)  and h4test > (h4dis-1000))):
    if h4test < h4dis:
        T4dis = T4dis + 0.1
    else:
        T4dis = T4dis -0.1
    HEOS.update(CoolProp.PT_INPUTS, P4dis, T4dis)
    h4test = HEOS.hmass()
T4dis = HEOS.T()
s4dis = HEOS.smass()

"""----Hot point lower stage----"""
h4hot = h4dis-((h4dis-h3)*compressor_loss)
P4hot = P4dis
h4test = 1000
T4hot = T4dis
while (not(h4test < (h4hot+1000)  and h4test > (h4hot-1000))):
    if h4test < h4hot:
        T4hot = T4hot + 0.1
    else:
        T4hot = T4hot - 0.1
    HEOS.update(CoolProp.PT_INPUTS, P4hot, T4hot)
    h4test = HEOS.hmass()
T4hot = HEOS.T()
s4hot = HEOS.smass()

"""---Condensation Point---"""
P4d = P4hot
HEOS.update(CoolProp.PQ_INPUTS, P4d, 1)
h4d = HEOS.hmass()

"""---Bubble point---"""
P4bub = P4d
HEOS.update(CoolProp.PQ_INPUTS, P4bub, 0)
h4bub = HEOS.hmass()
T4bub = HEOS.T()

"""----Subcooling point IHE Intermediate----"""
P4sub1 = P4bub - 1300
HEOS.specify_phase(CoolProp.iphase_liquid)
HEOS.update(CoolProp.PQ_INPUTS, P3, 0)
T3bub = HEOS.T()
T4sub1 = (T4bub-T3bub)*.3+T3bub

HEOS.update(CoolProp.PT_INPUTS, P3, T4sub1)
h4sub1 = HEOS.hmass()

"""----Subcooling point IHE Intermediate----"""
P4sub2 = P4sub1
h4sub2 = h4sub1 - (h2-h1)
h4test = 1000
T4sub2 = T4sub1
while (not(h4test < (h4sub2+1000)  and h4test > (h4sub2-1000))):
    if h4test < h4sub2:
        T4sub2 = T4sub2 + 0.1
    else:
        T4sub2 = T4sub2 - 0.1
    HEOS.specify_phase(CoolProp.iphase_liquid)
    HEOS.update(CoolProp.PT_INPUTS, P4sub2, T4sub2)
    h4test = HEOS.hmass()
T4sub2 = HEOS.T()
s4sub2 = HEOS.smass()

"""----Flash point----"""
h4f = h4sub2

"""---I point----"""
h5 = h4f
P5 = P1
HEOS.update(CoolProp.PQ_INPUTS, P5, 1)
HEOS.specify_phase(CoolProp.iphase_gas)
hg = HEOS.hmass()
HEOS.update(CoolProp.PQ_INPUTS, P5, 0)
HEOS.specify_phase(CoolProp.iphase_liquid)
hf = HEOS.hmass()
X5=(h5-hf)/(hg-hf) 

HEOS.specify_phase(CoolProp.iphase_twophase)
HEOS.update(CoolProp.PQ_INPUTS, P5, X5)
T5 = HEOS.T()
s5 = HEOS.smass()

"""---Intermediate expansion part---"""
P6 = P3
h6 = h4sub1
HEOS.update(CoolProp.PQ_INPUTS, P6, 1)
HEOS.specify_phase(CoolProp.iphase_gas)
hg6 = HEOS.hmass()
HEOS.update(CoolProp.PQ_INPUTS, P6, 0)
HEOS.specify_phase(CoolProp.iphase_liquid)
hf6 = HEOS.hmass()
X6=(h6-hf6)/(hg6-hf6)

HEOS.specify_phase(CoolProp.iphase_twophase)
HEOS.update(CoolProp.PQ_INPUTS, P6, X6)
T6 = HEOS.T()
s6 = HEOS.smass()

HEOS.build_phase_envelope("dummy")
PE = HEOS.get_phase_envelope_data()
PELabel = 'Propane, x = ' + str(propane_mass_fraction)
plt.plot([x / HEOS.molar_mass() for x in PE.hmolar_vap],PE.p, '-', label=PELabel)


h=[h1,h2,h3disintermediate,h3, h4dis, h4hot, h4d, h4bub, h4sub2, h5, h1,h2,h3disintermediate,h3, h4dis, h4hot, h4d, h4sub1, h6, h3]
P=[P1,P2,P3disintermediate,P3, P4dis, P4hot, P4d, P4bub, P4sub2, P5, P1,P2,P3disintermediate,P3, P4dis, P4hot, P4d, P4sub1, P6, P3]

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

plt.plot(h,P,'red')
plt.xlabel('Enthalpy [J/kg]')
plt.ylabel('Pressure [kPa]')
plt.yscale('log')
plt.ylim(5e4,10e6)
plt.xlim(200000,750000)
plt.title('Refrigerant cycle for R290/R600 Mixtures')
plt.legend(loc='lower right', shadow=True)
