from Simulation import*
    
def main(): 
    
    print(__name__)
    plt.figure(figsize=(20,5))
    plt.plot(data[0],label='Tair')
    plt.plot(data[1],label='Twall')
    plt.plot(SP_Sim,label='SP_Temperature')
    #plt.plot(T_outdoor_Sim,label='Toutdoor')
    plt.legend(loc='best')
    plt.show()

    #data = house(T_outdoor_Sim,Qinternal_Sim,Qsolar_Sim,SP_Sim,time_sim,CF,Rair_outdoor,Rair_wall,Cair,Cwall)    
# Using the special variable  
# __name__ 
if __name__=="__main__":
    main()
    	
    
