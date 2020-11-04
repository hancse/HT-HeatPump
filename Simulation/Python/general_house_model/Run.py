#import matplotlib.pyplot as plt
#import house as house
#import parameters as par
#import internal_heat_gain as ihg
#import Temperature_SP as sp
from house_model.house import *
from house_model.parameters import *
from house_model.internal_heat_gain import *
from house_model.Temperature_SP import *
from house_model.Total_Irrad import *
#from house_model import house
#import house_model

#Define Simulation time
days_Sim = 365                          #number of simulation days


time_sim      = time[0:days_Sim*24]
Qsolar_Sim    = Qsolar[0:days_Sim*24] 
Qinternal_Sim = Qinternal[0:days_Sim*24]
#Qinst_Sim     = Qinst_Sim[0:days_Sim*24][:,0]
T_outdoor_Sim = Toutdoor_p[0:days_Sim*24]

#Set point
SP_Sim=SP[0:days_Sim*24]
#CF=CF
#Rair_outdoor=Rair_outdoor
#Rair_wall=Rair_wall
#Cair=Cair
#Cwall=Cwall

data = house(T_outdoor_Sim,Qinternal_Sim,Qsolar,SP_Sim,time_sim,CF,Rair_outdoor,Rair_wall,Cair,Cwall)
