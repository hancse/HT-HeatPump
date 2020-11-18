# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 19:06:45 2020

@author: TrungNguyen
"""
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from weeks import read_week
import pandas as pd

def temp_sp_week_day(t1,t2,Night_T_SP,Day_T_SP,Wu_time,Work_time,back_home):
   
    """ Define Temperature SP for 1 days (24 hours).

    Args:
        t1:           Presence from [hour]
        t2:           Presence until [hour]
        Night_T_SP:   Set temperature of thermostat at night from time t2
        Day_T_SP:     Set wishes temperature of thermostat
        Wu_time:      Define wake up time in the morning, temperature set to 20
        duty_wu:
        Work_time:    Define time that people go to work.
        duty_w:
        back_home:    Define time that people back from work 18:00

    Returns:

    """

    # t1= 8                                #Presence from [hour]
    # t2= 23                               #Presence until [hour]

    # Night_T_SP=17                        # Set temperature of thermostat at night from time t2
    # Day_T_SP=20							 # Set wishes temperature of thermostat

    # Define Wake up time
    # Wu_time =7           				 # Define wake up time in the morning, temperature set to 20
    duty_wu = t2-Wu_time

    # Go to work time/ leave the house
    # Work_time = 8           			 # Define time that people go to work.
    duty_w   = t2-Work_time

    # Back to home
    # back_home = 18                       #Define time that people back from work 18:00
    duty_b   = t2-back_home

    # Creating profile

    days_hours   = 24                    #number_of_hour_in_oneday + start hour at 0
    days         = 365                   #number of simulation days
    periods      = 24*3600*days          #in seconds (day_periods*365 = years)
    pulse_width  = (t2-t1)/24            # % of the periods
    phase_delay  = t1                    #in seconds

    # temperature different between day and night.
    delta_T= Day_T_SP - Night_T_SP
    
	#-----------------------
    t= np.linspace(0,1,(days_hours*days)+1,endpoint=False)          #+1 start from 0 days=1
    temp1 = signal.square(2 * np.pi* days * t,duty=duty_wu/24)
    temp1 = np.clip(temp1, 0, 1)
	# add delay to array
    temp1=np.roll(temp1,Wu_time)

	#----------------
    t= np.linspace(0,1,(days_hours*days)+1,endpoint=False)          #+1 start from 0 days=1
    temp2 = signal.square(2 * np.pi* days * t,duty=duty_w/24)
    temp2 = np.clip(temp2, 0, 1)
	# add delay to array
    temp2=np.roll(temp2,Work_time)

	#___________
    t= np.linspace(0,1,(days_hours*days)+1,endpoint=False)          #+1 start from 0 days=1
    temp3 = signal.square(2 * np.pi* days * t,duty=duty_b/24)
    temp3 = np.clip(temp3, 0, 1)
	# add delay to array
    temp3=np.roll(temp3,back_home)

	# Calculate SP
    temp4=temp1-temp2+temp3
    SP_weekday=(temp4*delta_T)+Night_T_SP

    SP_weekday=SP_weekday[np.newaxis]
    SP_weekday=SP_weekday.T
    SP_weekday=np.delete(SP_weekday, -1, 0)
    
	#SP_weekday=SP_weekday.flatten()	
    
    return SP_weekday

def temp_sp_day_off(t1,t2,Night_T_SP,Day_T_SP,Wu_time,heating_off_time):
    
     # Creating profile

    days_hours   = 24                    #number_of_hour_in_oneday + start hour at 0
    days         = 365                   #number of simulation days
    periods      = 24*3600*days          #in seconds (day_periods*365 = years)
    pulse_width  = (t2-t1)/24            # % of the periods
    phase_delay  = t1                    #in seconds
    
    # temperature different between day and night.
    delta_T= Day_T_SP - Night_T_SP
		

    """Define wake up time """
    #W_time  = 10;           """wake up at 10: during day-off """
    duty_wo = t2-Wu_time        

    """Heating system turn off at 22:00"""
    #Sleep_time    = 22;           """Heating system turn off at 22:00 """
    duty_ho       = t2 - heating_off_time         


	#-----------------------
    t= np.linspace(0,1,(days_hours*days)+1,endpoint=False)          #+1 start from 0 days=1
    temp1 = signal.square(2 * np.pi* days * t,duty=duty_wo/24)
    temp1 = np.clip(temp1, 0, 1)
    # add delay to array
    temp1=np.roll(temp1,Wu_time)

	#----------------
    t= np.linspace(0,1,(days_hours*days)+1,endpoint=False)          #+1 start from 0 days=1
    temp2 = signal.square(2 * np.pi* days * t,duty=duty_ho/24)
    temp2 = np.clip(temp2, 0, 1)
    # add delay to array
    temp2=np.roll(temp2,heating_off_time)


	# Calculate SP

    temp4=temp1-temp2#+temp3
    SP_dayoff=(temp4*delta_T)+Night_T_SP

    SP_dayoff=SP_dayoff[np.newaxis]
    SP_dayoff=SP_dayoff.T
    SP_dayoff=np.delete(SP_dayoff, -1, 0)
	#SP_dayoff= SP_dayoff.flatten()
    return SP_dayoff

def SP_profile(df_week: pd.DataFrame,SP_weekday,SP_dayoff):
        
    temp2 = np.zeros((0,1))
    
    # A years profile loop
    for i in range(len(df_week)):
        
        temp  = SP_weekday[0:24]
        
        if df_week.Dag[i]=='zaterdag' or df_week.Dag[i] == 'zondag':
            #print(i)
            temp  = SP_dayoff[0:24]
            #print('------------------weekend------------------')
        for x in range (11):
            #print(x)
            if df_week.Datum[i] == df_week.holidays[x]:
                print(i)
                temp  = SP_dayoff[0:24]
        if i == 360:
            #plt.plot(SP_weekday[0:23])
            plt.plot(temp[7:24])

        temp2 = np.concatenate((temp2, temp))
 
    SP = temp2[0:8760].flatten()
    
    return SP


if __name__ == "__main__":
    
    df_weeks = read_week('NEN_data')
    week_day_setpoint = temp_sp_week_day(8, 23, 17, 20, 7, 8, 18)
    day_off_setpoint  = temp_sp_day_off(8, 23, 17, 20, 10, 22)
    print(df_weeks.Dag[364])
    #print(week_data.iloc[0]['Datum'] == week_data.iloc[0]['holidays'])
    SP =SP_profile(df_weeks,week_day_setpoint,day_off_setpoint)
    
    # Plot 48 hours
    plt.figure(figsize=(5, 5))
    #plt.plot(week_day_setpoint[0:48], label='setpoint')
    #plt.plot(day_off_setpoint[0:48], label='setpoint')
    plt.plot(SP[24*364+7:24*364+24], label='setpoint')
    #plt.ylabel('Temperature_SP (degC)')
    #plt.xlabel('time (h)')
    #plt.legend(loc='best')
    plt.show()