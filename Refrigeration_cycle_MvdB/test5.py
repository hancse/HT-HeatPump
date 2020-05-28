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
h = []
P = []
h1 = []
#for x0 in [0.02, 0.2, 0.4, 0.6, 0.8, 0.98]:
for P0 in range(130000, 750000, 100):
    h.append(PropsSI('H','T',273,'P',P0,'Propane[0.5]&IsoButane[0.5]'))
    P.append(P0)
    h1.append(PropsSI('H','T',273+10,'P',P0,'Propane[0.5]&IsoButane[0.5]'))
    
for x0 in [0.5]:
    HEOS.set_mass_fractions([x0, 1 - x0])
    try:
        HEOS.build_phase_envelope("dummy")
        PE = HEOS.get_phase_envelope_data()
        PELabel = 'Propane, x = ' + str(x0)
        plt.plot([x / HEOS.molar_mass() for x in PE.hmolar_vap], PE.p, '-', label=PELabel)
        
        plt.plot(h,P,'red')
        plt.plot(h1,P,'red')
        plt.xlabel('Enthalpy [J/kg]')
        plt.ylabel('Pressure [kPa]')
        plt.yscale('log')
        plt.ylim(10e4,10e6)
        plt.xlim(100000,750000)
        plt.title('Refrigerant cycle for Propane/IsoButane Mixtures')
        plt.legend(loc='lower right', shadow=True)
    except ValueError as VE:
        print(VE)    



    