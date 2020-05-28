import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import StateContainer

#setting the gasses
gas = 'R290[0.5]&R600[0.5]'
HEOS = CoolProp.AbstractState('HEOS', 'R290&R600')
HEOS.set_mass_fractions([0.5, 0.5])

#Setting the points
evaporation_temp = 0 + 273.15                
condensation_temp = 55 + 273.15  
superheat = 5                               
subcooling = 3.1  
isentropic_eff = 0.7
compressor_loss = 0.1


"""Evaporation"""

#PropsSi environment
T1 = evaporation_temp
P1 = PropsSI('P','T',T1,'Q',1,gas)
h1= PropsSI('H','T',T1,'Q',1,gas)    
s1= PropsSI('S','T',T1,'Q',1,gas)

#HEOS environment
HEOS.update(CoolProp.QT_INPUTS, 1, T1)
HEOS.specify_phase(CoolProp.iphase_gas)
HEOSP1 = HEOS.p()
HEOSh1 = HEOS.hmass()
HEOSs1 = HEOS.smass()

"""Superheating"""

#PropsSi environment
P2 = P1
T2 = superheat + evaporation_temp
h2= PropsSI('H','T|gas',T2,'P',P2,gas)    
s2= PropsSI('S','T|gas',T2,'P',P2,gas)

HEOS.update(CoolProp.PT_INPUTS, P1, T1)
HEOSP2 = HEOS.p()
HEOSh2 = HEOS.hmass()
HEOSs2 = HEOS.smass()

"""100% isentopic compression"""
#HEOS environment
T3d = condensation_temp
HEOS.update(CoolProp.QT_INPUTS, 1, T3d)
P3s = HEOS.p()
HEOSh3s = HEOS.hmass()
HEOSs3s = HEOS.smass()
T3s = condensation_temp
HEOS.update(CoolProp.PT_INPUTS, P3s, T3s)
s3s = HEOS.smass()

while (not(s3s < (s2+1)  and s3s > (s2-1))):
    if s3s < s2:
        T3s = T3s + 0.1
    else:
        T3s = T3s -0.1
    HEOS.update(CoolProp.PT_INPUTS, P3s, T3s)
    s3s = HEOS.smass()
h3s = HEOS.hmass()

"""Discharge point"""
h3dis = ((h3s-h2)/isentropic_eff) + h2
P3dis = P3s
HEOS.update(CoolProp.PT_INPUTS, P3dis, T3s)
T3dis = HEOS.T()
h3test = HEOS.hmass()

while (not(h3test < (h3dis+100)  and h3test > (h3dis-100))):
    if h3test < h3dis:
        T3dis = T3dis + 0.1
    else:
        T3dis = T3dis -0.1
    HEOS.update(CoolProp.PT_INPUTS, P3dis, T3dis)
    h3test = HEOS.hmass()
T3dis = HEOS.T()
s3dis = HEOS.smass()
#T3dis = PropsSI('T','H',h3dis,'P',P3dis,gas)

"""Hot point"""
h3hot = h3dis-((h3dis-h2)*compressor_loss)
P3hot = P3dis
h3test = 1000
T3hot = T3dis
while (not(h3test < (h3hot+1000)  and h3test > (h3hot-1000))):
    if h3test < h3hot:
        T3hot = T3hot + 0.1
    else:
        T3hot = T3hot - 0.1
    HEOS.update(CoolProp.PT_INPUTS, P3hot, T3hot)
    h3test = HEOS.hmass()
T3hot = HEOS.T()
s3hot = HEOS.smass()

"""condensation point"""
T3d = condensation_temp
P3d = P3dis
h3d = PropsSI('H','T|gas',T3d,'Q',1,gas)    
s3d = PropsSI('S','T|gas',T3d,'Q',1,gas)

"""Bubble point"""
P4b = P3dis
HEOS.update(CoolProp.PQ_INPUTS, P4b, 0)
HEOS.specify_phase(CoolProp.iphase_liquid)
T4b = HEOS.T()
h4b = HEOS.hmass()
h4b = PropsSI('H','T|liquid',T4b,'Q',0,gas)    
s4b = PropsSI('S','T|liquid',T3d,'Q',0,gas)

"""Subcooling point"""
P4sub = P4b
h4sub = h4b - (h2-h1) 
h4test = 1000
T4sub = T4b
while (not(h4test < (h4sub+1000)  and h4test > (h4sub-1000))):
    if h4test < h4sub:
        T4sub = T4sub + 0.1
    else:
        T4sub = T4sub - 0.1
    HEOS.update(CoolProp.PT_INPUTS, P4sub, T4sub)
    h4test = HEOS.hmass()
T4sub = HEOS.T()
s4sub = HEOS.smass()

"""Flash point"""
h4f = h4sub
#HEOS.update(CoolProp.HmassQ_INPUTS, h4f, 0)
#T4f = PropsSI('T','H|liquid',h4f,'Q',0,gas)

"""i point"""
h5 = h4f
P5 = P1
HEOS.update(CoolProp.PQ_INPUTS, P5, 1)
HEOS.specify_phase(CoolProp.iphase_gas)
hg = HEOS.hmass()
HEOS.update(CoolProp.PQ_INPUTS, P5, 0)
HEOS.specify_phase(CoolProp.iphase_liquid)
hf = HEOS.hmass()
X5=(h5-hf)/(hg-hf) 

HEOS.specify_phase(CoolProp.iphase_twophase)
HEOS.update(CoolProp.PQ_INPUTS, P5, X5)
T5 = HEOS.T()
s5 = HEOS.smass()

"""
T3d= condensation_temp                              #Bubble point
s3d= PropsSI('S','T',T3d,'Q',1,gas)
#T3s= PropsSI('T','P',P3s,'Q',1,gas)
P3d = PropsSI('P','T',T3d,'Q',1,gas)
#h3d = PropsSI('H','P',P3s,'S',s3s,gas) 
h3d = PropsSI('H','T',T3d,'Q',1,gas)   

HEOS.update(CoolProp.PT_INPUTS, P2, T2)
s2test = HEOS.smass()
kappa = HEOS.cpmolar()/HEOS.cvmolar()
T3stest =  ((P3d/P1)**((kappa-1)/kappa))*T2

while (not(s3stest < (s2+1)  and s3stest > (s2-1))):
    if s3stest < s2:
        T3stest = T3stest + 0.1
    else:
        T3stest = T3stest -0.1
    HEOS.update(CoolProp.PT_INPUTS, P3d, T3stest)
    s3stest = HEOS.smass()

h3stest = HEOS.hmass()
h3s = PropsSI('H','P',P3d,'S',s2,gas) 
s3s=  PropsSI('S','P',P3d,'H',h3s,gas)
"""
cycle_states = StateContainer()

cycle_states[" Evap" ,'H'] = h1
cycle_states[" Evap"]['S'] = s1
cycle_states[" Evap"][CoolProp.iP] = P1
cycle_states[" Evap",CoolProp.iT] = T1

cycle_states[" Sup " ,'H'] = h2
cycle_states[" Sup "]['S'] = s2
cycle_states[" Sup "][CoolProp.iP] = P2
cycle_states[" Sup ",CoolProp.iT] = T2

cycle_states[" Dis " ,'H'] = h3dis
cycle_states[" Dis "]['S'] = s3dis
cycle_states[" Dis "][CoolProp.iP] = P3dis
cycle_states[" Dis ",CoolProp.iT] = T3dis

cycle_states[" Hot " ,'H'] = h3hot
cycle_states[" Hot "]['S'] = "nan"
cycle_states[" Hot "][CoolProp.iP] = P3hot
cycle_states[" Hot ",CoolProp.iT] = "nan"

cycle_states[" Isen " ,'H'] = h3s
cycle_states[" Isen "]['S'] = s3s
cycle_states[" Isen "][CoolProp.iP] = P3s
cycle_states[" Isen ",CoolProp.iT] = T3s

cycle_states[" Cond " ,'H'] = h3d
cycle_states[" Cond "]['S'] = s3d
cycle_states[" Cond "][CoolProp.iP] = P3d
cycle_states[" Cond ",CoolProp.iT] = T3d

cycle_states[" Cond " ,'H'] = h3d
cycle_states[" Cond "]['S'] = s3d
cycle_states[" Cond "][CoolProp.iP] = P3d
cycle_states[" Cond ",CoolProp.iT] = T3d

cycle_states[" Bub " ,'H'] = h4b
cycle_states[" Bub "]['S'] = s4b
cycle_states[" Bub "][CoolProp.iP] = P4b
cycle_states[" Bub ",CoolProp.iT] = T4b

cycle_states[" Sub " ,'H'] = h4sub
cycle_states[" Sub "]['S'] = s4sub
cycle_states[" Sub "][CoolProp.iP] = P4sub
cycle_states[" Sub ",CoolProp.iT] = T4sub

cycle_states[" I " ,'H'] = h5
cycle_states[" I "]['S'] = s5
cycle_states[" I "][CoolProp.iP] = P5
cycle_states[" I ",CoolProp.iT] = T5

print(cycle_states)