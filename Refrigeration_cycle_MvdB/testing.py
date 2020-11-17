# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:23:49 2020

@author: Uitleen
"""
import CoolProp
from CoolProp.CoolProp import PropsSI
from prettytable import PrettyTable

t = PrettyTable(['Fraction Isobutane', 'Temperature'])

for i in range(0,11):
    frac = i/10
    t.add_row([frac, PropsSI('T','Q',1,'P',101325,'HEOS::R600a[%s]&R290[%s]'%(frac,1-frac))-273])

print(t)

