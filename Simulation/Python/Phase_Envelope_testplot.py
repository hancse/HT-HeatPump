import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots import StateContainer
import CoolProp
import numpy as np


CP.set_config_bool(CP.OVERWRITE_BINARY_INTERACTION, True)

#A = CP.apply_simple_mixing_rule('Propane', 'IsoButane', 'linear')
A = CP.apply_simple_mixing_rule('Propane', 'IsoButane', 'Lorentz-Berthelot')

#Evaporator chracteristics
gas = 'Propane&IsoButane'
isentropic_eff=0.7
evaporation_temp = 0 + 273.15                # R290
condensation_temp = 55 + 273.15              # Estimation from Aeronamic ??
superheat = 5                                # ambient temperature
subcooling = 3.1     


#D= CP.PropsSI('Dmass','T',300,'P',101325,'Propane[0.7]&IsoButane[0.3]')

#C= CP.PropsSI('C','T',300,'P',101325,'Propane[0.7]&IsoButane[0.3]')

#H= CP.PropsSI('H','T',300,'P',101325,'Propane[0.7]&IsoButane[0.3]')

HEOS = CP.AbstractState('HEOS', 'Propane&IsoButane')

#for x0 in [0.02, 0.2, 0.4, 0.6, 0.8, 0.98]:
for x0 in [0.7]:
    HEOS.set_mole_fractions([x0, 1 - x0])
    try:
        HEOS.build_phase_envelope("dummy")
        PE = HEOS.get_phase_envelope_data()
        PELabel = 'Propane, x = ' + str(x0)
        plt.plot(PE.T, PE.p, '-', label=PELabel)
    except ValueError as VE:
        print(VE)

#T= CP.PropsSI('T','H',28627,'P',101325,'Propane[0.7]&IsoButane[0.3]')
MapPT=[PE.T[0:-1],PE.p[0:-1]]        # create PT database
MapPTt=np.transpose(MapPT)


plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [Pa]')
plt.yscale('log')
plt.title('Phase Envelope for Propane/IsoButane Mixtures')
plt.legend(loc='lower right', shadow=True)
#plt.savefig('methane-ethane.pdf')
#plt.savefig('methane-ethane.png')
plt.show()

#T= CP.PropsSI('T','H',H,'P',101325,'Propane[0.7]&IsoButane[0.3]')
#PE.Temperature


P1=4.61950637e+05
T1 = evaporation_temp

# Create some test value
h1=PE.hmolar_vap[82]
h2=PE.hmolar_vap[100]

P2=PE.p[100]
h1=[h1,h2]
P1=[P1,P2]


plt.plot(PE.hmolar_vap,PE.p, '-', label=PELabel)
plt.plot(h1,P1,'red')
plt.show()
cycle_states = StateContainer()











