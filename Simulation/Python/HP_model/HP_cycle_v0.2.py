# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:42:21 2020

@author: Uitleen
"""
import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots import StateContainer


#Evaporator chracteristics
gas = 'R290'
isentropic_eff=0.7
evaporation_temp = 0 + 273.15
condensation_temp = 55 + 273.15
superheat = 5
subcooling = 3.1

#Compressor characteristics
expected_compressor_power = 4300
expected_massflow = 0.04506


#Aeronamic calculated values
Qdot_condensor_expected = 14.3e3
Expected_COP = 3.17


#after evaporator, not superheated
T1 = evaporation_temp
P1 = PropsSI('P','T',T1,'Q',1,gas)
h1= PropsSI('H','T',T1,'Q',1,gas)    
s1= PropsSI('S','T',T1,'Q',1,gas)

#Superheated Evap temp = 5C
T2 = superheat + evaporation_temp
P2 = P1
h2= PropsSI('H','T',T2,'Q',1,gas)    
s2= PropsSI('S','T',T2,'Q',1,gas)

#Calculating isentropic enthalpy
T_isen= condensation_temp
P_isen = PropsSI('P','T',condensation_temp,'Q',1,gas)
print(P_isen)
s_isen= s2
h_isen= PropsSI('H','S',s_isen,'P',P_isen,gas) 
print(h_isen)

#Temperature after correction for isentropic efficientcy
h_hot = (h_isen-h2)/0.7+h2
T_hot = PropsSI('T','H',h_hot,'P',P_isen,gas)
s_hot= PropsSI('S','T',T_hot,'Q',1,gas)
print(h_hot)
#temperature after correction for compressor efficiency
h_dis = ((h_isen-h2)*.1)+h_hot
print(h_dis)
T_dis = PropsSI('T','H',h_dis,'P',P_isen,gas)
s_dis= PropsSI('S','T',T_dis,'Q',1,gas)


P3 = P_isen
T3 = T_dis
h3= h_dis
s3=s_dis

#saturated Liquid line
T4 = condensation_temp
P4 = P_isen
s4 = PropsSI('S','T',T4,'Q',0,gas)
h4 = PropsSI('H','T',T4,'Q',0,gas) 



#Subcooling
T5 = condensation_temp - subcooling
P5 = PropsSI('P','T',T5,'Q',0,gas)
h5= PropsSI('H','T',T5,'Q',0,gas)    
s5= PropsSI('S','T',T5,'Q',0,gas)

T6 = evaporation_temp
h6=h5
P6 = P1
Q6 = PropsSI('Q','H',h6,'P',P6,gas)
h6= PropsSI('H','T',T6,'Q',Q6,gas)    
s6= PropsSI('S','T',T6,'Q',Q6,gas)

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

cycle_states[7,'H'] = h1
cycle_states[7]['S'] = s1
cycle_states[7][CoolProp.iP] = P1
cycle_states[7,CoolProp.iT] = T1

#Qdot_condensor_calculated = (h2-h4)*massflow

print(cycle_states)
print(Qdot_condensor_expected)
#print(Qdot_condensor_calculated)
print(Expected_COP)
#print(Qdot_condensor_calculated/compressor_power)


pp = PropertyPlot(gas, 'PH', unit_system='KSI',tp_limits='ACHP')
pp.calc_isolines(CoolProp.iQ, num=11)
pp.set_axis_limits([300,700,400,2500])
pp.draw_process(cycle_states)
pp.show()

pp = PropertyPlot(gas, 'TS', unit_system='KSI',tp_limits='ACHP')
pp.calc_isolines(CoolProp.iQ, num=15)
pp.set_axis_limits([1.2,2.6,250,360])
pp.draw_process(cycle_states)
pp.show()

pp = PropertyPlot(gas, 'PH', unit_system='KSI',tp_limits='ACHP')
pp.calc_isolines(CoolProp.iQ, num=11)
pp.set_axis_limits([300,700,400,2500])
pp.draw_process(cycle_states)
pp.show()
#plot = PropertyPlot('R290', 'ph')
#plot.calc_isolines()
#plot.show()