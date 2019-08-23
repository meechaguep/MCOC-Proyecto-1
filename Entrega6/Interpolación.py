# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 22:57:56 2019

@author: Vicente
"""
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

Temp=open('Temperatura.txt')
Time=open('Tiempo.txt')
temperatura=[]
tiempo=[]
for i in Temp:
    a=float(i)
    temperatura.append(a)
for j in Time:
    b=float(j)
    tiempo.append(b)
f = interp1d(tiempo, temperatura)
xnwe=np.linspace(0,10)
plt.plot(tiempo, f(tiempo), '-')
