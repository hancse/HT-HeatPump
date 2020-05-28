# -*- coding: utf-8 -*-
"""
Created on Thu May 7 15:27:45 2020

@author: Maarten van den Berg
"""
import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots import StateContainer


#Evaporator chracteristics
gas = 'R290'
isentropic_eff=0.7
evaporation_temp = 0 + 273.15                # R290
condensation_temp = 55 + 273.15              # Estimation from Aeronamic
superheat = 5
superheat_economizer = 5                                # ambient temperature
subcooling_evap = 5                          # estimation from Aeronamic (how to know or calculate the value)
subcooling_economizer = 5
intermediate_temp = 30 + 273.15
subcooling = 10

#Compressor characteristics
compressor_power = 4300                      #
massflow = 0.04506                           #


#Aeronamic calculated values
Qdot_condensor_expected = 14.3e3
Expected_COP = 3.17


#after evaporator, not superheated
T1 = evaporation_temp
P1 = PropsSI('P','T',T1,'Q',1,gas)
h1= PropsSI('H','T',T1,'Q',1,gas)    
s1= PropsSI('S','T',T1,'Q',1,gas)

T2 = superheat + evaporation_temp
P2 = P1
#h2= PropsSI('H','T',T2,'Q',1,gas)    
#s2= PropsSI('S','T',T2,'Q',1,gas)
h2= PropsSI('H','T|gas',T2,'P',P2,gas)    
s2= PropsSI('S','T|gas',T2,'P',P2,gas)

#Isentropic enthalpy base on P,S (Min work)
P3s=  PropsSI('P','T',intermediate_temp,'Q',1,gas)
h3s = PropsSI('H','P',P3s,'S',s2,gas) 
T3s=  PropsSI('T','P',P3s,'H',h3s,gas) 
s3s=  PropsSI('S','P',P3s,'H',h3s,gas)

#isentropic_eff=0.7
h3 = ((h3s-h2)/isentropic_eff) + h2
P3 = P3s
T3 = PropsSI('T','P',P3,'H',h3,gas)
#h3 = PropsSI('H','T|gas',T3,'P',P3,gas)    
s3 = PropsSI('S','T',T3,'P',P3,gas)

#intermdediate stage
P8 = P3s
T8 = intermediate_temp + superheat_economizer
h8 = PropsSI('H','T|gas',T8,'P',P8,gas)
s8 = PropsSI('S','T',T8,'P',P8,gas)

#intermdediate stage
P4 = P3s
h4 = ((h3s-h8)/2)+h8
T4 = PropsSI('T','P',P4,'H',h4,gas)
s4 = PropsSI('S','T',T8,'P',P8,gas)

#Condensor level
P5 = P5s =  PropsSI('P','T',condensation_temp,'Q',1,gas)
h5s = PropsSI('H','P',P5s,'S',s8,gas)
T5s=  PropsSI('T','P',P5s,'H',h5s,gas) 
s5s=  PropsSI('S','P',P5s,'H',h5s,gas)

#isentropic_eff=0.7
h5 = ((h5s-h4)/isentropic_eff) + h4
T5 = PropsSI('T','P',P5,'H',h5,gas)
#h3 = PropsSI('H','T|gas',T3,'P',P3,gas)    
s5 = PropsSI('S','T',T5,'P',P5,gas)

#Bubble Point
P6 = P5
T6 = condensation_temp
h6 = PropsSI('H','T',T6,'Q',0,gas)
#s6 = PropsSI('S','T',T6,'P',P6,gas)

#Bubble Point
P6 = P5
T6 = condensation_temp
h6 = PropsSI('H','T',T6,'Q',0,gas)
#s6 = PropsSI('S','T',T6,'P',P6,gas)

#Subcooling
P6s=P5
T6s = condensation_temp - subcooling
s6s= PropsSI('S','T',T6s,'P',P6s,gas)
h6s = PropsSI('H','T',T6s,'P',P6s,gas) 
#PropsSI('T','P',P2s,'Q',1,gas)

T8 = 273.15 + 35
P8 = P4
h8 = h6
                           #Pressure drop evap (evap =P1 , P5=P1+ pressure drop)
#s5= PropsSI('S','T',T5,'Q',X5,gas)




cycle_states = StateContainer()

cycle_states[1,'H'] = h1
cycle_states[1]['S'] = s1
cycle_states[1][CoolProp.iP] = P1
cycle_states[1,CoolProp.iT] = T1

cycle_states[2,'H'] = h2
cycle_states[2]['S'] = s2
cycle_states[2][CoolProp.iP] = P2
cycle_states[2,CoolProp.iT] = T2

cycle_states[3,'H'] = h3
cycle_states[3]['S'] = s3
cycle_states[3][CoolProp.iP] = P3
cycle_states[3,CoolProp.iT] = T3

cycle_states[4,'H'] = h4
cycle_states[4]['S'] = s4
cycle_states[4][CoolProp.iP] = P4
cycle_states[4,CoolProp.iT] = T4

cycle_states[5,'H'] = h5
cycle_states[5]['S'] = s5
cycle_states[5][CoolProp.iP] = P5
cycle_states[5,CoolProp.iT] = T5

cycle_states[6,'H'] = h6
cycle_states[6]['S'] = s6
cycle_states[6][CoolProp.iP] = P6
cycle_states[6,CoolProp.iT] = T6

cycle_states[7,'H'] = h6s
cycle_states[7]['S'] = s6s
cycle_states[7][CoolProp.iP] = P6s
cycle_states[7,CoolProp.iT] = T6s

cycle_states[8,'H'] = h7
#cycle_states[8]['S'] = s7
cycle_states[8][CoolProp.iP] = P7
cycle_states[8,CoolProp.iT] = T7

cycle_states[10,'H'] = h1
cycle_states[10]['S'] = s1
cycle_states[10][CoolProp.iP] = P1
cycle_states[10,CoolProp.iT] = T1

cycle_states[11,'H'] = h2
cycle_states[11]['S'] = s2
cycle_states[11][CoolProp.iP] = P2
cycle_states[11,CoolProp.iT] = T2

cycle_states[12,'H'] = h3
cycle_states[12]['S'] = s3
cycle_states[12][CoolProp.iP] = P3
cycle_states[12,CoolProp.iT] = T3

cycle_states[13,'H'] = h4
cycle_states[13]['S'] = s4
cycle_states[13][CoolProp.iP] = P4
cycle_states[13,CoolProp.iT] = T4

cycle_states[14,'H'] = h5
cycle_states[14]['S'] = s5
cycle_states[14][CoolProp.iP] = P5
cycle_states[14,CoolProp.iT] = T5

cycle_states[15,'H'] = h6
cycle_states[15]['S'] = s6
cycle_states[15][CoolProp.iP] = P6
cycle_states[15,CoolProp.iT] = T6

cycle_states[16,'H'] = h8
#cycle_states[16]['S'] = s6
cycle_states[16][CoolProp.iP] = P8
cycle_states[16,CoolProp.iT] = T8

cycle_states[17,'H'] = h4
#cycle_states[17]['S'] = s6
cycle_states[17][CoolProp.iP] = P4
cycle_states[17,CoolProp.iT] = T4
"""
Qdot_condensor_calculated = (h3-h4)*massflow
Calculated_COP=Qdot_condensor_calculated/compressor_power
print(cycle_states)
print("Qdot_condensor_expected:",Qdot_condensor_expected)
print("Qdot_condensor_calculated:",Qdot_condensor_calculated)
print("Expected_COP:",Expected_COP)
print("Calculated_COP:",Calculated_COP)
"""



pp = PropertyPlot(gas, 'PH', unit_system='KSI',tp_limits='ACHP')
pp.calc_isolines(CoolProp.iQ, num=11)
pp.set_axis_limits([300,700,400,2500])
pp.draw_process(cycle_states)
pp.show()
