# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 12:05:19 2020

@author: TrungNguyen
"""
# Defining main function
import house 				as hs
import parameters			as par
import internal_heat_gain 	as hg
import Temperature_SP		as sp
import Total_Irrad			as irrad
import matplotlib.pyplot 	as plt

def main():
	
	# Define Simulation time 
	days_Sim = 20                          #number of simulation days
	time_sim      = par.time[0:days_Sim*24]
	Qsolar_Sim    = par.Qsolar[0:days_Sim*24]
	#Qsolar_Sim    = Qsolar[0:days_Sim*24]*0 
	 
	Qinternal_Sim = hg.Qinternal[0:days_Sim*24]
	#Qinst_Sim     = Qinst_Sim[0:days_Sim*24][:,0]
	T_outdoor_Sim = par.Toutdoor[0:days_Sim*24]

	#Set point

	SP_Sim		 = sp.SP[0:days_Sim*24]
	CF			 = par.CF
	Rair_outdoor = par.Rair_outdoor
	Rair_wall	 = par.Rair_wall
	Cair		 = par.Cair
	Cwall		 = par.Cwall

	data = hs.house(T_outdoor_Sim,Qinternal_Sim,Qsolar_Sim,SP_Sim,time_sim,CF,Rair_outdoor,Rair_wall,Cair,Cwall)
    
	#______plot the results________________
	plt.figure(figsize=(20,5))
	plt.plot(data[0],label='Tair')
	plt.plot(data[1],label='Twall')
	plt.plot(SP_Sim,label='SP_Temperature')
	#plt.plot(T_outdoor_Sim,label='Toutdoor')
	plt.legend(loc='best')
	plt.show()
	

main()
	
#main()
