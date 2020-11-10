# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 12:05:19 2020

@author: TrungNguyen
"""
# Defining main function
from house import *
from parameters import *
from internal_heat_gain import *
from Temperature_SP import *
from Total_Irrad import *

# Define Simulation time 
days_Sim = 20                          #number of simulation days
time_sim      = time[0:days_Sim*24]
Qsolar_Sim    = Qsolar[0:days_Sim*24]
#Qsolar_Sim    = Qsolar[0:days_Sim*24]*0 
 
Qinternal_Sim = Qinternal[0:days_Sim*24]
#Qinst_Sim     = Qinst_Sim[0:days_Sim*24][:,0]
T_outdoor_Sim = Toutdoor_p[0:days_Sim*24]

#Set point

SP_Sim=SP[0:days_Sim*24]
CF=CF
Rair_outdoor=Rair_outdoor
Rair_wall=Rair_wall
Cair=Cair
Cwall=Cwall

data = house(T_outdoor_Sim,Qinternal_Sim,Qsolar_Sim,SP_Sim,time_sim,CF,Rair_outdoor,Rair_wall,Cair,Cwall)    
