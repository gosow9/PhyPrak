# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 11:33:43 2021

@author: fritz
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from gauss import error





print('------------------------------------')
print('Experiment 1')
print('------------------------------------')

V_1 = np.array([200, 300, 400, 500, 600, 700, 800, 900, 980]) #start with 200, +100
t_1 = np.ones_like(V_1)*30#so 1% is reached
N_1 = np.array([0, 0, 1826, 2379, 2535, 2452, 2501, 2402, 2516])
n_1 = N_1/t_1

del_V = np.ones_like(V_1)*10
del_t = 0
del_N = np.sqrt(N_1)

def get_steigung(x):
    return np.polyfit(x[0,3:], x[1, 3:], 1)[0]

del_pol = error(get_steigung, np.array([V_1, N_1]), np.array([del_V, del_N]))

poly = np.polyfit(V_1[3:], N_1[3:], 1)
lin = lambda V: poly[0]*V + poly[1]
V_poly = np.linspace(500, 980)


plt.figure(figsize=(6.4, 4))
plt.plot(V_poly, lin(V_poly), label = 'Geiger plateau fit')
plt.errorbar(V_1, N_1, del_N, del_V, '-x', label = 'Measurements')
plt.ylabel('Number of decays, N')
plt.xlabel('Voltage, in Volts')
plt.legend()
plt.savefig('plateau.PNG')
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
del_N_2 = np.sqrt(N_2)
N_bg = 50
del_N_bg = np.sqrt(N_bg)
Neff_2 = N_2-N_bg
del_Neff_2 = np.sqrt(del_N_2**2 + del_N_bg**2)

d = 101 #mm
r = 18.05 #mm Radius Stahlblende #pm 0.025
del_r = 0.5*0.025
del_d = 0.5

get_theta = lambda x: np.arctan(x[0]/2/x[1]) #x = (r, d)
theta = get_theta(np.array([r, d])) #deg Halbwinkel
del_theta = error(get_theta, np.array([r, d]), np.array([del_r, del_d]))


get_eps = lambda theta: np.sin(theta/2)**2
eps = get_eps(theta)
del_eps = error(get_eps, np.array([theta]), np.array([del_theta]))

get_aktiv = lambda x: x[0]/x[1]/x[2] #x = (N, t, eps)
aktiv_2 = get_aktiv(np.array([Neff_2, t_2, eps]))
del_aktiv_2 = error(get_aktiv, np.array([Neff_2, t_2, eps]), np.array([del_Neff_2, 0., del_eps[0]]))

get_aktiv_tilde = lambda aktiv: aktiv/2*(1/0.45 + 1/0.9)

aktiv_tilde = get_aktiv_tilde(aktiv_2) #Bq
del_aktiv_tilde = error(get_aktiv_tilde, np.array([aktiv_2]), np.array([del_aktiv_2]))[0]
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

print('N = ', N_2, ' +- ', np.sqrt(N_2))
print('N_bg = ', N_bg, ' +- ', np.sqrt(N_bg))
print('N_eff = ', Neff_2, ' +- ', del_Neff_2)
print('Halbwinkel: ', np.degrees(theta), '+-', np.degrees(del_theta),  'deg')
print('Akzeptanz: ', eps, '+-', del_eps)
print('Aktivität: ', aktiv_2, '+-', del_aktiv_2, 'Bq')
print('Aktivität ohne Folie: ', aktiv_tilde, '+-', del_aktiv_tilde, 'Bq')


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
rho = 2.69 #mg cm^-3
x = np.array([0, 0.3, 0.5, 1.1, 1.6, 2.,  2.3, 2.5, 2.8, 3., 3.5, 4.])*10 #cm
E_verlust = 125e-3


rho_x = rho*x #g cm^-2
Ex = np.array([0, 0.07, 0.09, 0.15, 0.2, 0.23, 0.26, 0.28, 0.29, 0.3, 0.32, 0.34]) #abgelesen abb. 4
E_max_1 = 0.32

t_3 = 90
N_3 = np.array([7390, 4606, 3124, 1496, 689, 308, 205, 128, 82, 67, 51, 53])
Neff_3 = N_3 - N_bg
aktiv_3 = Neff_3/t_3/eps
E = np.array([])
Neff_n = Neff_3**(1/n)
x_max = 3.5
rho_max = rho*x_max

E_max_q = E_max_1 + E_verlust
mu =0.
E_max_mu = 0.


plt.figure(figsize=(6.4, 4))
plt.plot(x, N_3)
plt.plot(x, np.ones_like(x)*N_bg)
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

print('E_max abgelesen: ', E_max_1, 'MeV')
print('E_max_quelle: ', E_max_q, 'MeV')
















