# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 11:33:43 2021

@author: fritz
"""

import numpy as np
import matplotlib.pyplot as plt


print('------------------------------------')
print('Experiment 1')
print('------------------------------------')

V = np.array([1.]) #start with 200, +100
t = np.array([1.]) #so 1% is reached
N = np.array([1.])

plt.figure(figsize=(12, 8))
plt.plot(N, V, '-')
plt.xlabel('Number of decays, N')
plt.ylabel('Voltage, V')
plt.show()

i = 0
V_ideal = V[i]
print('V_ideal: ', V_ideal)
print('del_N/del_V: ', V_ideal/N[i])


print('\n------------------------------------')
print('Experiment 2')
print('------------------------------------')
t = 0
N = 0
N_bg = 0
Neff = N-N_bg

d = 101 #mm
r = 0 #mm Radius Stahlblende
theta = 0 #deg Halbwinkel
eps = np.sin(theta/2)**2

aktiv = 0 #Bq
aktiv_tilde = 0 #Bq

print('Aktivität verwendete Probe: ',aktiv)
print('Aktivität ohne Folie: ', aktiv_tilde)

print('\n------------------------------------')
print('Experiment 3')
print('------------------------------------')

















