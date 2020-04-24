# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:27:45 2020

@author: TrungNguyen
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
superheat = 5                                # ambient temperature
subcooling = 3.1                             # estimation from Aeronamic (how to know or calculate the value)

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

#Dew point calculation
T3d= condensation_temp                              #Bubble point
s3d= PropsSI('S','T',T3d,'Q',1,gas)
#T3s= PropsSI('T','P',P3s,'Q',1,gas)
P3d = PropsSI('P','T',T3d,'Q',1,gas)
#h3d = PropsSI('H','P',P3s,'S',s3s,gas) 
h3d = PropsSI('H','T',T3d,'Q',1,gas)                #Dew_Enthalpy : check graph

#Isentropic enthalpy base on P,S (Min work)
P3s=P3d
h3s = PropsSI('H','P',P3d,'S',s2,gas) 
T3s=  PropsSI('T','P',P3d,'H',h3s,gas) 
s3s=  PropsSI('S','P',P3d,'H',h3s,gas)

# Hot enthalpy and hot gas calculation

#hot gas and hot enthalpy temp ??? ask

hhot=(h3s-h2)
Phot=P3s
Thot= 72.8 +273.15     # Temp?????
hhot = PropsSI('H','T',Thot,'P',P3d,gas)
shot = PropsSI('S','T',Thot,'P',P3d,gas) 


#isentropic_eff=0.7
h3= ((h3s-h2)/isentropic_eff) + h2
P3 = P3d
T3 = PropsSI('T','P',P3,'H',h3,gas)
#h3 = PropsSI('H','T|gas',T3,'P',P3,gas)    
s3 = PropsSI('S','T',T3,'P',P3,gas)


#isentropic_e= (h3s-h2)/(h3-h2)
#print(isentropic_e)

#saturated Liquid line bub

P4 = P3
P4=P4 - 12*1000                             #Pressure drop assumption
T4 =PropsSI('T','P',P4,'Q',0,gas)
#T4= T3d
#P4= PropsSI('P','T',T4,'Q',0,gas)
h4= PropsSI('H','P',P4,'Q',0,gas)    
s4= PropsSI('S','P',P4,'Q',0,gas)

#Subcooling
P4s=P4
T4s = condensation_temp - subcooling
s4s= PropsSI('S','T',T4s,'P',P4s,gas)
h4s = PropsSI('H','T',T4s,'P',P4s,gas) 
#PropsSI('T','P',P2s,'Q',1,gas) 

#flash point
h4f=h4s
T4f = T4s
s4f= PropsSI('S','T',T4f,'Q',0,gas)
P4f = PropsSI('P','T',T4f,'Q',0,gas) 

#PropsSI('T','P',P2s,'Q',1,gas) 

# i point
T5 = T1
h5=h4s
P5 = P1
P5= P5 + 8*1000                              #Pressure drop evap (evap =P1 , P5=P1+ pressure drop)
hg=PropsSI('H','T',T1,'Q',1,gas)
hf=PropsSI('H','T',T1,'Q',0,gas)
X5=(h5-hf)/(hg-hf)                           # X =Quality
s5= PropsSI('S','T',T5,'Q',X5,gas)           # X in vapor


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

cycle_states[4,'H'] = hhot
cycle_states[4]['S'] = shot
cycle_states[4][CoolProp.iP] = P3
cycle_states[4,CoolProp.iT] = Thot

cycle_states[5,'H'] = h3s
cycle_states[5]['S'] = s3s
cycle_states[5][CoolProp.iP] = P3
cycle_states[5,CoolProp.iT] = T3s

cycle_states[6,'H'] = h3d
cycle_states[6]['S'] = s3d
cycle_states[6][CoolProp.iP] = P3
cycle_states[6,CoolProp.iT] = T3d

cycle_states[7,'H'] = h4
cycle_states[7]['S'] = s4
cycle_states[7][CoolProp.iP] = P4
cycle_states[7,CoolProp.iT] = T4

cycle_states[8,'H'] = h4s
cycle_states[8]['S'] = s4s
cycle_states[8][CoolProp.iP] = P4s
cycle_states[8,CoolProp.iT] = T4s

cycle_states[9,'H'] = h4f
cycle_states[9]['S'] = s4f
cycle_states[9][CoolProp.iP] = P4f
cycle_states[9,CoolProp.iT] = T4f

cycle_states[10,'H'] = h5
cycle_states[10]['S'] = s5
cycle_states[10][CoolProp.iP] = P5
cycle_states[10,CoolProp.iT] = T5


cycle_states[11,'H'] = h1
cycle_states[11]['S'] = s1
cycle_states[11][CoolProp.iP] = P1
cycle_states[11,CoolProp.iT] = T1



Qdot_condensor_calculated = (h3-h4)*massflow
Calculated_COP=Qdot_condensor_calculated/compressor_power
print(cycle_states)
print("Qdot_condensor_expected:",Qdot_condensor_expected)
print("Qdot_condensor_calculated:",Qdot_condensor_calculated)
print("Expected_COP:",Expected_COP)
print("Calculated_COP:",Calculated_COP)


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