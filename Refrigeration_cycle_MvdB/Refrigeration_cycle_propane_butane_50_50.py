import CoolProp
from CoolProp.Plots import StateContainer
import matplotlib.pyplot as plt

#setting the gasses
#propane_mass_fraction = 0.5
HEOS = CoolProp.AbstractState('HEOS', 'IsoButane')
HEOS.build_phase_envelope("dummy")
#HEOS.set_mass_fractions([0.5, 0.5])

#Setting the points
evaporation_temp = 0 + 273.15                
condensation_temp = 75 + 273.15  
superheat = 5                               
isentropic_eff = 0.70
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


"""----100% isentopic compression----"""
T3d = condensation_temp
HEOS.update(CoolProp.QT_INPUTS, 1, T3d)
P3s = HEOS.p()

#create inital guesses for T3s and s3s
T3s = condensation_temp
HEOS.update(CoolProp.PT_INPUTS, P3s, T3s)
s3s = HEOS.smass()

#Using a loop to find real values for H3s and T3s
while (not(s3s < (s2+1)  and s3s > (s2-1))):
    if s3s < s2:
        T3s = T3s + 0.1
    else:
        T3s = T3s -0.1
    HEOS.update(CoolProp.PT_INPUTS, P3s, T3s)
    s3s = HEOS.smass()
h3s = HEOS.hmass()

"""----Discharge point----"""
h3dis = ((h3s-h2)/isentropic_eff) + h2
P3dis = P3s
HEOS.update(CoolProp.PT_INPUTS, P3dis, T3s)
T3dis = HEOS.T()
h3test = HEOS.hmass()

while (not(h3test < (h3dis+1000)  and h3test > (h3dis-1000))):
    if h3test < h3dis:
        T3dis = T3dis + 0.1
    else:
        T3dis = T3dis -0.1
    HEOS.update(CoolProp.PT_INPUTS, P3dis, T3dis)
    h3test = HEOS.hmass()
T3dis = HEOS.T()
s3dis = HEOS.smass()

"""----Hot point----"""
h3hot = h3dis-((h3dis-h2)*compressor_loss)
P3hot = P3dis
h3test = 1000
T3hot = T3dis
while (not(h3test < (h3hot+1000)  and h3test > (h3hot-1000))):
    if h3test < h3hot:
        T3hot = T3hot + 0.1
    else:
        T3hot = T3hot - 0.1
    HEOS.update(CoolProp.PT_INPUTS, P3hot, T3hot)
    h3test = HEOS.hmass()
T3hot = HEOS.T()
s3hot = HEOS.smass()

"""----condensation point----"""
T3d = condensation_temp
P3d = P3dis
HEOS.update(CoolProp.QT_INPUTS, 1, T3d)
h3d = HEOS.hmass()   
s3d = HEOS.smass()

"""----Bubble point----"""
P4b = P3dis
HEOS.update(CoolProp.PQ_INPUTS, P4b, 0)
HEOS.specify_phase(CoolProp.iphase_liquid)
T4b = HEOS.T()
h4b = HEOS.hmass()
s4b = HEOS.smass()

"""----Subcooling point----"""
P4sub = P4b - 1300
h4sub = h4b - (h2-h1) 
h4test = 1000
T4sub = T4b
while (not(h4test < (h4sub+1000)  and h4test > (h4sub-1000))):
    if h4test < h4sub:
        T4sub = T4sub + 0.1
    else:
        T4sub = T4sub - 0.1
    HEOS.specify_phase(CoolProp.iphase_liquid)
    HEOS.update(CoolProp.PT_INPUTS, P4sub, T4sub)
    h4test = HEOS.hmass()
T4sub = HEOS.T()
s4sub = HEOS.smass()


"""----Flash point----"""
h4f = h4sub


"""----i point----"""
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

"""----------------------------------------------------"""

cycle_states = StateContainer()

cycle_states[" Evap" ,'H'] = h1
cycle_states[" Evap"]['S'] = s1
cycle_states[" Evap"][CoolProp.iP] = P1
cycle_states[" Evap",CoolProp.iT] = T1

cycle_states[" Sup " ,'H'] = h2
cycle_states[" Sup "]['S'] = s2
cycle_states[" Sup "][CoolProp.iP] = P2
cycle_states[" Sup ",CoolProp.iT] = T2

cycle_states[" Dis " ,'H'] = h3dis
cycle_states[" Dis "]['S'] = s3dis
cycle_states[" Dis "][CoolProp.iP] = P3dis
cycle_states[" Dis ",CoolProp.iT] = T3dis

cycle_states[" Hot " ,'H'] = h3hot
cycle_states[" Hot "]['S'] = s3hot
cycle_states[" Hot "][CoolProp.iP] = P3hot
cycle_states[" Hot ",CoolProp.iT] = T3hot

cycle_states[" Isen " ,'H'] = h3s
cycle_states[" Isen "]['S'] = s3s
cycle_states[" Isen "][CoolProp.iP] = P3s
cycle_states[" Isen ",CoolProp.iT] = T3s

cycle_states[" Cond " ,'H'] = h3d
cycle_states[" Cond "]['S'] = s3d
cycle_states[" Cond "][CoolProp.iP] = P3d
cycle_states[" Cond ",CoolProp.iT] = T3d

cycle_states[" Cond " ,'H'] = h3d
cycle_states[" Cond "]['S'] = s3d
cycle_states[" Cond "][CoolProp.iP] = P3d
cycle_states[" Cond ",CoolProp.iT] = T3d

cycle_states[" Bub " ,'H'] = h4b
cycle_states[" Bub "]['S'] = s4b
cycle_states[" Bub "][CoolProp.iP] = P4b
cycle_states[" Bub ",CoolProp.iT] = T4b

cycle_states[" Sub " ,'H'] = h4sub
cycle_states[" Sub "]['S'] = s4sub
cycle_states[" Sub "][CoolProp.iP] = P4sub
cycle_states[" Sub ",CoolProp.iT] = T4sub

cycle_states[" I " ,'H'] = h5
cycle_states[" I "]['S'] = s5
cycle_states[" I "][CoolProp.iP] = P5
cycle_states[" I ",CoolProp.iT] = T5

print(cycle_states)


PE = HEOS.get_phase_envelope_data()
#PELabel = 'Propane, x = ' + str(propane_mass_fraction)
plt.plot([x / HEOS.molar_mass() for x in PE.hmolar_vap],PE.p, '-')


h=[h1,h2,h3dis, h4sub, h5, h1]
P=[P1,P2,P3dis, P4sub, P5, P1]


plt.plot(h,P,'red')
plt.xlabel('Enthalpy [J/kg]')
plt.ylabel('Pressure [Pa]')
plt.yscale('log')
plt.ylim(10e4,10e6)
plt.xlim(250000,750000)
plt.title('Refrigerant cycle for R290/R600 Mixtures')
plt.legend(loc='lower right', shadow=True)

Massflow_calculated = condensor_capacity/(h3hot-h4b)
Calculated_COP=(h3hot-h4b)/(h3dis-h2)
compressor_power = (h3dis-h2)*Massflow_calculated
print("Massflow:",Massflow_calculated)
print("Calculated_COP:",Calculated_COP)
print("Calculated_Compressor power:",compressor_power)