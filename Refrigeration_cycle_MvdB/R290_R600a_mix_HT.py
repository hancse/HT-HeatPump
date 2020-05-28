import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI


#Evaporator chracteristics
isentropic_eff=0.7
evaporation_temp = 0 + 273.15                # R290
condensation_temp = 55 + 273.15              # Estimation from Aeronamic ??
superheat = 5                                # ambient temperature
subcooling = 10     


#D= CP.PropsSI('Dmass','T',300,'P',101325,'Propane[0.7]&IsoButane[0.3]')

#C= CP.PropsSI('C','T',300,'P',101325,'Propane[0.7]&IsoButane[0.3]')

#H= CP.PropsSI('H','T',300,'P',101325,'Propane[0.7]&IsoButane[0.3]')

HEOS = CP.AbstractState('HEOS', 'Propane&IsoButane')

#for x0 in [0.02, 0.2, 0.4, 0.6, 0.8, 0.98]:
for x0 in [0.99]:
    HEOS.set_mass_fractions([x0, 1 - x0])

    HEOS.build_phase_envelope("dummy")
    PE = HEOS.get_phase_envelope_data()
    PELabel = 'Propane, x = ' + str(x0)
    plt.plot([x / HEOS.molar_mass() for x in PE.hmolar_vap], PE.p, '-', label=PELabel)
    
    T1 = evaporation_temp
    P1 = PropsSI('P','T',T1,'Q',1,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    h1= PropsSI('H','T',T1,'Q',1,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    s1= PropsSI('S','T',T1,'Q',1,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    
    T2 = superheat + evaporation_temp
    P2 = P1
    #h2= PropsSI('H','T',T2,'Q',1,gas)    
    #s2= PropsSI('S','T',T2,'Q',1,gas)
    h2= PropsSI('H','T|gas',T2,'P',P2,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    s2= PropsSI('S','T|gas',T2,'P',P2,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    
    #Dew point calculation
    T3d= condensation_temp                              #Bubble point
    s3d= PropsSI('S','T',T3d,'Q',1,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    #T3s= PropsSI('T','P',P3s,'Q',1,gas)
    P3d = PropsSI('P','T',T3d,'Q',1,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    #h3d = PropsSI('H','P',P3s,'S',s3s,gas) 
    h3d = PropsSI('H','T',T3d,'Q',1,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))   

    #Isentropic enthalpy base on P,S (Min work)
    P3s=P3d
    s3s = 2000
    kappa = HEOS.cpmolar()/HEOS.cvmolar()
    T3s =  ((P3d/P1)**((kappa-1)/kappa))*T2
    
    while (not(s3s < (s2+1)  and s3s > (s2-1))):
        if s3s < s2:
            T3s = T3s + 0.1
        else:
            T3s = T3s -0.1
        HEOS.update(CP.PT_INPUTS, P3d, T3s)
        s3s = HEOS.smass()

    
    h3s = HEOS.hmass()
    h3= ((h3s-h2)/isentropic_eff) + h2
    P3 = P3d
    #T3 = PropsSI('T','P',P3,'H',h3,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    #h3 = PropsSI('H','T|gas',T3,'P',P3,gas)    
    #s3 = PropsSI('S','T',T3,'P',P3,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    
    P4 = P3d
    P4=P4 - 12*1000                             #Pressure drop assumption
    T4 =PropsSI('T','P',P4,'Q',0,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    #T4= T3d
    #P4= PropsSI('P','T',T4,'Q',0,gas)
    h4= PropsSI('H','P',P4,'Q',0,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))    
    s4= PropsSI('S','P',P4,'Q',0,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    
    #Subcooling
    P4s=P4
    T4s = condensation_temp - subcooling
    s4s= PropsSI('S','T',T4s,'P',P4s,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    h4s = PropsSI('H','T',T4s,'P',P4s,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0)) 
    #PropsSI('T','P',P2s,'Q',1,gas) 
    
    #flash point
    h4f=h4s
    T4f = T4s
    s4f= PropsSI('S','T',T4f,'Q',0,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    P4f = PropsSI('P','T',T4f,'Q',0,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    
    T5 = T1
    h5=h4s
    P5 = P1
    P5= P5 + 8*1000                              #Pressure drop evap (evap =P1 , P5=P1+ pressure drop)
    hg=PropsSI('H','T',T1,'Q',1,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    hf=PropsSI('H','T',T1,'Q',0,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))
    X5=(h5-hf)/(hg-hf)                           # X =Quality
    s5= PropsSI('S','T',T5,'Q',X5,'Propane[%s]&IsoButane[%s]'%(x0, 1 - x0))           # X in vapor
    
    h=[h1,h2,h3s, h3, h4, h4s, h4f, h5, h1]
    P=[P1,P2,P3s, P3, P4, P4s, P4f, P5, P1]
    
    plt.plot(h,P,'red')
    plt.xlabel('Enthalpy [J/kg]')
    plt.ylabel('Pressure [kPa]')
    plt.yscale('log')
    plt.ylim(10e4,10e6)
    plt.xlim(250000,750000)
    plt.title('Refrigerant cycle for Propane/IsoButane Mixtures')
    plt.legend(loc='lower right', shadow=True)


plt.show()
#plt.savefig('methane-ethane.pdf')
#plt.savefig('methane-ethane.png')

#T= CP.PropsSI('T','H',H,'P',101325,'Propane[0.7]&IsoButane[0.3]')
#PE.Temperature
















