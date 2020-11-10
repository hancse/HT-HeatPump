'''house model base on 2R2C model

'''

from scipy.integrate import odeint
import numpy as np # linear algebra


def house(T_outdoor,Q_internal,Q_solar,SP_T,time_sim,CF,Rair_outdoor,Rair_wall,Cair,Cwall):
    """Compute air and wall tempearature inside the house.

        Args:
            T_outdoor (array) :  Outdoor tempearature in degree C
            Q_internal (array):  Internal heat gain in w.

            Q_solar (array)   :  Solar irradiation on window [W]
            SP_T (array)	  :  Setpoint tempearature from thermostat.
            Qinst (array)	  :  instant heat from heat source such as HP or boiler [W].

            time_sim (array)  :  simulation time
			
			CF				  : 
			Rair_outdoor	  :
			Rair_wall	      :
			Cair			  :
			Cwall             :

        Returns:
            tuple :  Tuple containing (Tair, Twall):

                - Tair (float):   air temperature inside the house in degree C.
                - Twall (float):  wall temperature inside the house in degree C
    """

		#Define model
    def model(x,t,T_outdoor,Q_internal,Q_solar,SP_T,Qinst,CF,Rair_outdoor,Rair_wall,Cair,Cwall):
        '''sub function making a model for odeint:
	
	       Args:
	 
	       x,t: ode input function 
	         func : callable(x, t, ...) or callable(t, x, ...)
			 Computes the derivative of y at t.
	         If the signature is ``callable(t, y, ...)``, then the argument
	         `tfirst` must be set ``True``.
		
	
    '''

	# States :
		
        Tair            = x[0]   # Temperature Buffer Tank (K)
        Twall           = x[1]   # Return Temperature to Floor (K)
		#integal         = x[2]
		#Qinstdt         = x[2]
		
	#     err      = SP_T-Tair
	#     integaldt= err
	#     integald= np.clip(integaldt, 0, 5)
	#     Qinst    = kP*(err) + ki*integal # kd*Tair.dt() # PID form
	#     Qinst=np.clip(Qinst, 0, 7000)
	#     #m.Equation(Integl.dt()==err )

		#___________
        Tairdt   = ((T_outdoor-Tair)/Rair_outdoor + (Twall-Tair)/Rair_wall  + Qinst + Q_internal + CF*Q_solar)/Cair
        Twalldt  = ((Tair-Twall)/Rair_wall + (1-CF)*Q_solar)/Cwall
		
		# Return array.astype(np.float64).ravel()
        return [Tairdt,Twalldt]
	# Initial Conditions for the States

    Tair0    = 20   
    Twall0   = 20

    y0 = [Tair0,Twall0]

    t = time_sim           # Define Simulation time with sampling time
    #print(len(t))

    Tair  = np.ones(len(t))*Tair0
    Twall  = np.ones(len(t))*Twall0
    consumption =  np.ones(len(t))
    kp= 7000  
    for i in range(len(t)-1):
        
        err=SP_T[i+1] - Tair[i]
        Qinst=err*kp
        Qinst=np.clip(Qinst, 0, 7000)
		
        inputs = (T_outdoor[i],Q_internal[i],Q_solar[i],SP_T[i],Qinst,CF,Rair_outdoor,Rair_wall,Cair,Cwall)
        #print(i)
        ts = [t[i],t[i+1]]
        y = odeint(model,y0,ts,args=inputs)

        Tair[i+1]  = y[-1][0]
        Twall[i+1] = y[-1][1]
		#integal[i+1]         = y[-1][2]

		# Adjust initial condition for next loop
        
        y0 = y[-1]
        consumption[i] = Qinst
    return Tair, Twall, consumption
#if __name__ == '__main__':
#	print('test')
#	house()