import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import StateContainer
import matplotlib.pyplot as plt

#setting the gasses
propane_mass_fraction = 0.5
HEOS = CoolProp.AbstractState('HEOS', 'Propane&R600')
HEOS.set_mass_fractions([0.5, 0.5])

#Setting the points
evaporation_temp = 0 + 273.15                
condensation_temp = 55 + 273.15
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
P3sintermediate = ((P4-P1)/2)+P1

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

HEOS.build_phase_envelope("dummy")
PE = HEOS.get_phase_envelope_data()
PELabel = 'Propane, x = ' + str(propane_mass_fraction)
plt.plot([x / HEOS.molar_mass() for x in PE.hmolar_vap],PE.p, '-', label=PELabel)


h=[h1,h2,h3disintermediate,h3, h4dis, h4hot]
P=[P1,P2,P3disintermediate,P3, P4dis, P4hot]


plt.plot(h,P,'red')
plt.xlabel('Enthalpy [J/kg]')
plt.ylabel('Pressure [kPa]')
plt.yscale('log')
plt.ylim(10e4,10e6)
plt.xlim(250000,750000)
plt.title('Refrigerant cycle for R290/R600 Mixtures')
plt.legend(loc='lower right', shadow=True)
