# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 16:02:55 2020

@author: Uitleen
"""

import CoolProp
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots import SimpleCompressionCycle
pp = PropertyPlot('R290', 'PH', unit_system='EUR',tp_limits='ACHP')
pp.calc_isolines(CoolProp.iQ, num=11)
cycle = SimpleCompressionCycle('R290', 'PH', unit_system='EUR')
T0 = 273
pp.state.update(CoolProp.QT_INPUTS,0.0,T0-10)
p0 = pp.state.keyed_output(CoolProp.iP)
T2 = 273+55
pp.state.update(CoolProp.QT_INPUTS,1.0,T2+15)
p2 = pp.state.keyed_output(CoolProp.iP)
pp.calc_isolines(CoolProp.iT, [T0-273.15,T2-273.15], num=2)
cycle.simple_solve(T0, p0, T2, p2, 0.76, SI=True)
cycle.steps = 50
sc = cycle.get_state_changes()
print(sc)
pp.draw_process(sc)
import matplotlib.pyplot as plt
plt.close(cycle.figure)
pp.show()