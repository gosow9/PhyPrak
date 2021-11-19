# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 11:33:43 2021

@author: fritz
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime




print('------------------------------------')
print('Experiment 1')
print('------------------------------------')

V_1 = np.array([200, 300, 400, 500, 600, 700, 800, 900, 980]) #start with 200, +100
t_1 = np.ones_like(V_1)*30#so 1% is reached
N_1 = np.array([0, 0, 1826, 2379, 2535, 2452, 2501, 2402, 2516])
poly = np.polyfit(V_1[3:], N_1[3:], 1)
lin = lambda V: poly[0]*V + poly[1]
V_poly = np.linspace(500, 980)


plt.figure(figsize=(12, 8))
plt.plot(V_poly, lin(V_poly))
plt.errorbar(V_1, N_1, None, 10, '-o')
plt.ylabel('Number of decays, N')
plt.xlabel('Voltage, V')
plt.show()

i = 0
V_ideal = 740
N_ideal = 2467
spann = N_ideal/V_ideal
print('V_ideal: ', V_ideal, 'V')
print('Spannungsabhängikeit: ', poly[0])


print('\n------------------------------------')
print('Experiment 2')
print('------------------------------------')
t_2 = 90
N_2 = 7453
N_bg = 50
Neff_2 = N_2-N_bg

d = 101 #mm
r = 18.05 #mm Radius Stahlblende #pm 0.025
theta = (np.tan(r/2/d)) #deg Halbwinkel

eps = np.sin(theta/2)**2

aktiv_2 = Neff_2/t_2/eps #Bq
aktiv_tilde = aktiv_2/2*(100/45+100/90) #Bq
t_half_sr = 28.8
t_half_y = 64
tau_sr = t_half_sr/np.log(2)
tau_y = t_half_y/np.log(2)
tau = tau_sr+tau_y/24/365 #years

then = datetime(1950, 1, 1)
now = datetime(2021, 11, 19)
diff = now-then
diff = diff.total_seconds()/60/60/24/365

A0 = aktiv_tilde/np.exp(-diff/tau)
t_10k = -tau*np.log(10e3/aktiv_tilde)


print('Halbwinkel: ', np.degrees(theta), 'deg')
print('Akzeptanz: ', eps)
print('Aktivität: ', aktiv_2, 'Bq')
print('Aktivität ohne Folie: ', aktiv_tilde, 'Bq')


print('\n------------------------------------')
print('Experiment 3')
print('------------------------------------')

D = 5.5e-3
d_SR = 2.8e-8
d_Y = 2.7e-9


t_dosis = D/(aktiv_tilde*(d_SR + d_Y))
print('1/e-Zeit Sr: ', tau_sr, 'y')
print('1/e-Zeit Y: ', tau_y, 'h')
print('Aktivität 1950: ', A0, 'Bq')
print('10k Bq in: ', t_10k, 'y')
print('Dosis erreicht nach: ', t_dosis, 's')


print('\n------------------------------------')
print('Experiment 4')
print('------------------------------------')
n = 1
rho = 2.69 #g cm^-3
x = np.array([0, 0.3, 0.5, 1.1, 1.6, 2.,  2.3, 2.5, 2.8, 3., 3.5, 4.])
rho_x = rho*x
t_3 = 90*np.ones_like(x)
N_3 = np.array([7390, 4606, 3124, 1496, 689, 308, 205, 128, 82, 67, 51, 53])
Neff_3 = N_3 - N_bg
aktiv = np.array([])
E = np.array([])
Neff_n = Neff_3**(1/n)
x_max = 0.
E_max =0.
E_verlust = 0.
E_max_q = 0.
mu =0.
E_max_mu = 0.


plt.figure(figsize=(6.4, 4))
plt.plot(x, N_3)
plt.yscale('log')
plt.xlabel('Thickness x of aluminium in cm')
plt.ylabel('Number of counts in 90s')
plt.show()

plt.figure(figsize=(6.4, 4))
plt.plot(x, N_3-N_bg)
plt.yscale('log')
plt.xlabel('Thickness x of aluminium in cm')
plt.ylabel('Number of counts in 90s')
plt.show()



print('E_max: ', E_max)















