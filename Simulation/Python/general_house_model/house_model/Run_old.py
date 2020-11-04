import matplotlib.pyplot as plt
#import house as house
import parameters as par
import internal_heat_gain as ihg
import Temperature_SP as sp
from house_model.house import house

#Define Simulation time
days_Sim = 25                          #number of simulation days

time_sim      = par.time[0:days_Sim*24]
Qsolar_Sim    = par.Qsolar[0:days_Sim*24] 
Qinternal_Sim = ihg.Qinternal[0:days_Sim*24]
#Qinst_Sim     = Qinst_Sim[0:days_Sim*24][:,0]
T_outdoor_Sim = par.Toutdoor_p[0:days_Sim*24]

#Set point
SP_Sim=sp.SP[0:days_Sim*24]

data = house.house(T_outdoor_Sim,Qinternal_Sim,Qsolar_Sim,SP_Sim,0,time_sim)