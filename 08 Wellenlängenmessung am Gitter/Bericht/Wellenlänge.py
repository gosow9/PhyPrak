# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 16:57:52 2021

@author: fritz
"""

import numpy as np
from gauss import error

print('------------------------------------')
print('Experiment 1')
print('------------------------------------')


g = 1/15000*25.4 #1/mm

b = np.array([230, 270, 280, 300, 305, 360, 420, 450])
del_b = np.ones_like(b)*7.5

ref1 = np.array([382, 447, 471, 492, 501, 587, 656, 706])

a = np.ones_like(b)*1000 #mm
del_a = np.ones_like(b)*5



get_tg_phi = lambda x: x[1]/x[0] #x = (a, b)
tg_phi = get_tg_phi(np.array([a, b]))
del_tg_phi = error(get_tg_phi, np.array([a, b]), np.array([del_a, del_b]))


get_phi = lambda tg: np.arctan(tg)
phi1 = np.arctan(tg_phi)
del_phi1 = error(get_phi, tg_phi, del_tg_phi)


get_lam1 = lambda phi: g*np.sin(phi) #x=(g, phi)
lam1 = get_lam1(phi1)
del_lam1 = error(get_lam1, phi1, del_phi1)

print('Wellenlänge Lambda = ', lam1*10e5, '+-\n', del_lam1*10e5, 'nm')


print('\n------------------------------------')
print('Experiment 2')
print('------------------------------------')

sec = lambda x: x/60


A1_1 = np.deg2rad(np.array([167+sec(20), 166+sec(15), 164+sec(20), 164+sec(11), 162+sec(24), 161+sec(17), 161+sec(12), 159+sec(36), 157+sec(10)]))
A2_1 = np.deg2rad(np.array([195+sec(24), 196+sec(33), 198+sec(26), 198+sec(35), 200+sec(27), 201+sec(29), 201+sec(33), 203+sec(10), 205+sec(37)]))

del_A = np.deg2rad(np.ones_like(A1_1)*sec(0.5))

A1_2 = np.deg2rad(np.array([153+sec(9), 150+sec(44), 146+sec(18), 145+sec(55), 141+sec(41), 138+sec(57), 138+sec(45), 134+sec(33), 127+sec(29)]))
A2_2 = np.deg2rad(np.array([210+sec(12), 212+sec(35), 217+sec(9), 217+sec(31), 221+sec(51), 224+sec(40), 224+sec(53), 229+sec(13), 236+sec(39)]))





get_phi2 = lambda A: 0.5*(A[0]-A[1])
phi2_1 = get_phi2(np.array([A2_1, A1_1]))
del_phi2_1 = error(get_phi2, np.array([A2_1, A1_1]), np.array([del_A, del_A]))

phi2_2 = get_phi2(np.array([A2_2, A1_2]))
del_phi2_2 = error(get_phi2, np.array([A2_2, A1_2]), np.array([del_A, del_A]))

lam2_1 = get_lam1(phi2_1)
del_lam2_1 = error(get_lam1, phi2_1, del_phi2_1)

get_lam2 = lambda phi: g*np.sin(phi)/2
lam2_2 = get_lam2(phi2_2)
del_lam2_2 = error(get_lam2, phi2_2, del_phi2_2)


print('1. Ordnung:')
print('Lambda = ', lam2_1*10e5, '+-\n', del_lam2_1*10e5, 'nm')

print('\n2. Ordnung:')
print('Lambda = ', lam2_2*10e5, '+-\n', del_lam2_2*10e5, 'nm')









