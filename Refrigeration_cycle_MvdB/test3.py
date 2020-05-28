# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:27:45 2020

@author: TrungNguyen
"""
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots import StateContainer


#Evaporator chracteristics
HEOS = CP.AbstractState('HEOS', 'Propane&IsoButane')
HEOS.set_mole_fractions([0.3, 0.7])

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
P1 = PropsSI('P','T',T1,'Q',1,'Propane&IsoButane')
h1= PropsSI('H','T',T1,'Q',1,'Propane&IsoButane')    
s1= PropsSI('S','T',T1,'Q',1,'Propane&IsoButane')