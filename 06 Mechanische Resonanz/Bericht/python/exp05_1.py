# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 12:24:54 2021

@author: Cedric
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq, fsolve, curve_fit
import matplotlib
from gauss import error
# matplotlib.use("pgf")
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })



"""
Experiment 1 bestimmung von Dämpfungskonstanten
"""
print("**********************************************************")
print("Experiment 1")
print("----------------------------------------------------------\n")

A0 = 110
T1 = 47.41 / 24
T1_err = 0.4/24
T2 = 27.60 / 14
T2_err = 0.4/14
T3 = 15.97 / 8
T3_err = 0.4/8
omega1 = 2*np.pi/(T1)
omega2 = 2*np.pi/(T2)
omega3 = 2*np.pi/(T3)
print(omega1, omega2, omega3)
PT1 = np.array([0, 3, 6, 9, 12, 15, 18, 21, 24])*T1  # Number of period
A11 = np.array([110, 84, 59, 42, 27, 19, 11, 7, 3])
A12 = np.array([110, 84, 60, 40, 29, 19, 12, 6, 2])
A1 = (A11+A12)/2   # Amplitude
lnA1 = np.log(A1/A1[0])
PT2 = np.array([0, 1, 2, 4, 6, 8, 10, 12, 14])*T2  # Number of periods
A2 = np.array([110, 97, 80, 56, 38, 25, 15, 10, 6])  # Amplitude
lnA2 = np.log(A2/A2[0])
PT3 = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])*T3  # Number of periods
A3 = np.array([110, 85, 64, 48, 36, 26, 19, 14, 10])  # Amplitude
lnA3 = np.log(A3/A3[0])

xerr = np.ones(A1.shape)*0.4 # +- 0,4 seconds time error
yerr = np.ones(A1.shape)*2   # +- 2 Degrees error
yerr[0] = 0.5
PT1x = np.linspace(0, PT1[-1],100)
PT2x = np.linspace(0, PT2[-1],100)
PT3x = np.linspace(0, PT3[-1],100)

alpha1, pcov1 = curve_fit(lambda t, a: 110*np.exp(-a*t),  PT1,  A1, 0.0)
expon1 = 110*np.exp(-alpha1[0]*PT1x)


alpha2 = curve_fit(lambda t, a: 110*np.exp(-a*t),  PT2,  A2, 0.0)[0]
expon2 = 110*np.exp(-alpha2[0]*PT2x)

alpha3 = curve_fit(lambda t, a: 110*np.exp(-a*t),  PT3,  A3, 0.0)[0]
expon3 = 110*np.exp(-alpha3[0]*PT3x)

print("alpha1=", alpha1)
print("alpha2=", alpha2)
print("alpha3=", alpha3)

plt.figure(figsize=(6.4, 4))
plt.ylabel(r'Amplitude A, $[A] = \deg$')
plt.xlabel(r'Time $[s] = \sec$')
plt.errorbar(PT1, A1, yerr=yerr, fmt=".", color ="red", capsize=3, label=r'Measured amplitudes dampening $I_1$')
plt.plot(PT1x, expon1, "--", color="red", label=r'Modell $A_0 e^{-\alpha_1 t}$ with $\alpha_1$ '+"= {:.2}".format(alpha1[0]))
plt.errorbar(PT2, A2, yerr=yerr, fmt=".", color ="blue", capsize=3, label=r'Measured amplitudes dampening $I_2$')
plt.plot(PT2x, expon2, "--", color="blue", label=r'Modell $A_0 e^{-\alpha_2 t}$ with $\alpha_2$ '+"= {:.2}".format(alpha2[0]))
plt.errorbar(PT3, A3, yerr=yerr, fmt=".", color ="green", capsize=3, label=r'Measured amplitudes dampening $I_3$')
plt.plot(PT3x, expon3, "--", color="green", label=r'Modell $A_0 e^{-\alpha_3 t}$ with $\alpha_3$ '+"= {:.3}".format(alpha3[0]))
plt.legend()
plt.savefig("damp1.PNG", dpi=600)


err_1 = np.ones(9)*2/A0

plt.figure(figsize=(6.4, 4))
plt.yscale("log")
plt.plot(PT1, A1/A0, ".",color="red", label="Dämpfung 1")
plt.errorbar(PT1, A1/A0, yerr=err_1, fmt=".", color="red", capsize=3, label=r'Measured amplitudes')
plt.plot(PT1x, expon1/A0, label="Model 1 ")
plt.plot(PT2, A2/A0, ".", label="Dämpfung 2")
plt.plot(PT2x, expon2/A0, label="Model 2")
plt.plot(PT3, A3/A0, ".", label="Dämpfung 3")
plt.plot(PT3x, expon3/A0, label="Model 13")


"""
Experiment 2 Eichung des Tachometers
"""
print("\n**********************************************************")
print("Experiment 2")
print("----------------------------------------------------------\n")
def get_C(x):
    return 2*np.pi*x[0]/(x[1]*x[2])

N_eich = np.array([30, 30, 70, 50])
del_N = np.zeros_like(N_eich)

t_eich = np.array([66.50, 51.59, 119.52, 108.59])
del_t = np.ones_like(t_eich)*0.2

V_tacho = np.array([1.157, 1.5, 1.5, 1.179])
del_V = np.ones_like(V_tacho)*0.0005

om_eich = 2*np.pi*N_eich/t_eich
C_eich = om_eich/V_tacho
C = np.mean(C_eich)
C_sx = error(get_C, np.array([N_eich, t_eich, V_tacho]), np.array([del_N, del_t, del_V]))
del_C = np.mean(C_sx)
print('C =', C, ' +- ', del_C, ' 1/Vs')





"""
Experiment 3 Erzwungene Schwingung (Resonanzkurven)
"""
print("\n**********************************************************")
print("Experiment 3")
print("----------------------------------------------------------\n")

om_test = np.linspace(0.4, 0.7)


get_omega = lambda V: V*C


#Dämpfung 1
V_t_1 = np.array([1.001, 1.037, 1.076, 1.111, 1.150, 1.189, 1.225, 1.265, 1.300, 1.341, 1.375, 1.412, 1.451, 1.489]) #in Volts
A_1 = np.array([7, 9, 9, 10, 15, 22, 32, 52, 83, 32, 17, 12, 10, 9])
om_1 = get_omega(V_t_1)
A_1_max = max(A_1)

poly_11 = np.polyfit(get_omega(np.array([1.265, 1.3])), np.array([52, 83]), 1)
poly_12 = np.polyfit(get_omega(np.array([1.3, 1.341])), np.array([83, 32]), 1)
z11 = (-poly_11[1]+A_1_max/np.sqrt(2))/poly_11[0]
z12 = (-poly_12[1]+ A_1_max/np.sqrt(2))/poly_12[0]
alpha_1=(abs(z12-z11)/2)

#Dämpfung 2
V_t_2 = np.array([1.002, 1.038, 1.075, 1.112, 1.151, 1.187, 1.225, 1.265, 1.301, 1.340, 1.375, 1.413, 1.450, 1.488])
A_2 = np.array([7, 8, 8, 12, 14, 20, 28, 44, 49, 27, 18, 14, 10, 9])
om_2 = get_omega(V_t_2)
A_2_max = max(A_2)

poly_21 = np.polyfit(get_omega(np.array([1.225, 1.265])), np.array([28, 44]), 1)
poly_22 = np.polyfit(get_omega(np.array([1.301, 1.340])), np.array([49, 27]), 1)
z21 = (-poly_21[1]+A_2_max/np.sqrt(2))/poly_21[0]
z22 = (-poly_22[1]+ A_2_max/np.sqrt(2))/poly_22[0]
alpha_2=(abs(z22-z21)/2)


#Dämpfung 3
V_t_3 = np.array([1.082, 1.130, 1.150, 1.187, 1.225, 1.263, 1.299, 1.338, 1.375, 1.413, 1.451, 1.487])
A_3 = np.array([9, 10, 12, 16, 23, 29, 27, 21, 17, 12, 10, 9])
om_3 = get_omega(V_t_3)
A_3_max = max(A_3)

poly_31 = np.polyfit(get_omega(np.array([1.187, 1.225])), np.array([16, 23]), 1)
poly_32 = np.polyfit(get_omega(np.array([1.338, 1.413])), np.array([21, 12]), 1)
z31 = (-poly_31[1]+A_3_max/np.sqrt(2))/poly_31[0]
z32 = (-poly_32[1]+ A_3_max/np.sqrt(2))/poly_32[0]
alpha_3=(abs(z32-z31)/2)

delta_V = 0.0005
get_del_omega = lambda C, del_C, V, del_V: np.sqrt((V*del_C)**2 + (C*del_V)**2)
del_om_1 = get_del_omega(C, del_C, V_t_1, delta_V)
del_om_2 = get_del_omega(C, del_C, V_t_2, delta_V)
del_om_3 = get_del_omega(C, del_C, V_t_3, delta_V)


plt.figure(figsize=(12, 8))
plt.errorbar(om_1, A_1, np.ones_like(om_1)*2, None, '--o', label = 'Dampening 1', capsize=3)
#plt.plot(np.array([z11, z12]), np.array([A_1_max, A_1_max])/np.sqrt(2), '-b')

plt.errorbar(om_2, A_2, np.ones_like(om_2)*2, None, '--o', label = 'Dampening 2', capsize=3)
#plt.plot(np.array([z21, z22]), np.array([A_2_max, A_2_max])/np.sqrt(2), '-g')

plt.errorbar(om_3, A_3, np.ones_like(om_3)*2, None, '--o', label = 'Dampening 3', capsize=3)
#plt.plot(np.array([z31, z32]), np.array([A_3_max, A_3_max])/np.sqrt(2), '-r')

plt.xlabel('Angular frequency $\omega$, $[\omega] = s^{-1}$')
plt.ylabel('Amplitude A, $[A] = \deg$')
plt.legend()
plt.savefig("resonance.PNG")


print('alpha 1: ', alpha_1)
print('alpha 2: ', alpha_2)
print('alpha 3: ', alpha_3)

print('delta alpha 1: ', abs(alpha_1-alpha1))
print('delta alpha 2: ', abs(alpha_2-alpha2))
print('delta alpha 3: ', abs(alpha_3-alpha3))

plt.show()


# T1m = T1.mean()
# T2m = T2.mean()
# dT1 = 0.4 / 5
# dT2 = 0.4 / 3
# print("Mean T1 = {}, std of T1 = {}".format(T1m, T1.std()))
# print("Mean T2 = {}, std of T2 = {}".format(T2m, T2.std()))
# Ms = 3.431             # Disc weight in kg +- 0.5 gramm
# Rs = 0.249/2       # Disc radius in m
# dM = 0.0005
# dR = 0.0005
# Os = moment_disc(Ms, Rs)
# D = torsion_disc(Os, T1.mean(), T2.mean())
# dOs = err_moment_disc(Ms, Rs, dM, dR, 1)
# dD = err_torsion_disc(Os, T1m, T2m, dOs, dT1, dT2, 5)
# print(r'Phi disc = {} +- {}'.format(Os, dOs))
# print(r'D disc = {} +- {}'.format(D, dD))
# print("**********************************************************\n\n")


# """
# Experiment 2 Satz von Steiner
# """
# print("Experiment 2")
# print("----------------------------------------------------------\n")


# Ta_5 = np.array([20.34, 20.38, 20.40])/3
# Ta_10 = np.array([22.50, 22.38, 22.50])/3
# Ta_15 = np.array([25.72, 25.69, 25.65])/3
# Ta_20 = np.array([29.71, 29.75, 29.69])/3
# Ta_25 = np.array([34.15, 33.97, 34.22])/3
# Ta = np.array([T1.mean(), Ta_5.mean(), Ta_10.mean(),
#                Ta_15.mean(), Ta_20.mean(), Ta_25.mean()])
# Mq = np.array([0.611, 0.611, 0.611, 0.605, 0.605, 0.605])  # Quad weight in kg
# d = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25])
# M = Mq.mean()
# c = 0.6
# b = 0.01174
# dc = 0.5/1000
# dd = dc*2
# db = 0.2/10000
# dTa = 0.4/3

# Oe = moment_ta(Ta, D)             # Experimental moment of inertia
# Obar = moment_quad(M, b, c)  # Calculated momentof inertia
# Oc = steiner(Obar, d, M)

# dOe = err_moment_ta(Ta, D, dTa, dD, 3)
# dOs = err_moment_quad(M, b, c, dM, db, dc)
# dOc = err_steiner(d, M, dOs, dd, dM)

# p = np.polyfit(d, Ta, 6)
# d_err = np.linspace(0, 0.25, 100)
# poly = np.poly1d(p)
# Ta_err = poly(d_err)
# Oe_err = moment_ta(Ta_err, D)

# dOe_err= err_moment_ta(Ta_err, D, dTa, dD, 3)

# # ----------------------------------------------------------------------------

# print("Mean 1 bar = {}kg, std of bars = {}, Two bars = {}kg".format(
#     Mq.mean(), Mq.std(), 2*Mq.mean()))
# print("Theta experimental = {} +- {}".format(Oe, dOe))
# print("Theta calculated = {} +- {}".format(Oc, dOc))
# print("**********************************************************\n\n")
# plt.figure()
# plt.plot(d, Ta, ".")
# plt.plot(d_err, Ta_err)



# plt.figure(figsize=(6.4, 4))
# plt.fill_between(d_err**2, Oe_err - dOe_err, Oe_err + dOe_err, color='red', alpha=0.1, label=r'Approximated errorband of $\theta_a$')
# plt.plot(d**2, Oc, label=r'Calculated moment of inertia $\theta_c$')
# #plt.plot(d_err**2, Oe_err, label=r'Calculated $\Theta_c$')
# #plt.plot(d**2, Oe, "or", label=r'Measured $\Theta_m$')
# plt.errorbar(d**2, Oe, xerr=None, yerr=dOe, fmt=".", color ="red", capsize=3, label=r'Measured moment of inertia $\theta_a$')
# plt.legend()
# plt.ylabel(r'Moment of inertia $\theta$')
# plt.xlabel(r'Squared displacement of bars $d^2$ in $m^2$')
# plt.savefig("steiner.pgf")


# """
# Experiment 3 Ellipse
# """
# print("Experiment 3")
# print("----------------------------------------------------------\n")
# Ts_0 = np.array([17.28, 17.37, 17.34])/5
# Ts_10 = np.array([17.65, 17.63, 17.66])/5
# Ts_20 = np.array([18.25, 18.28, 18.31])/5
# Ts_30 = np.array([19.34, 19.31, 19.25])/5
# Ts_40 = np.array([20.44, 20.40, 20.47])/5
# Ts_50 = np.array([21.59, 21.59, 21.53])/5
# Ts_60 = np.array([22.57, 22.46, 22.50])/5
# Ts_70 = np.array([23.31, 23.31, 23.35])/5
# Ts_80 = np.array([23.97, 23.85, 23.84])/5
# Ts_90 = np.array([23.97, 24.03, 24.06])/5
# # 0.2 Grad error
# phi = np.radians(np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90]))
# print(np.cos(phi)**2)

# Ts_std = np.array([Ts_0.std(), Ts_10.std(), Ts_20.std(), Ts_30.std(),
#                    Ts_40.std(), Ts_50.std(), Ts_60.std(), Ts_70.std(),
#                    Ts_80.std(), Ts_90.std()])
# Ts = np.array([Ts_0.mean(), Ts_10.mean(), Ts_20.mean(), Ts_30.mean(),
#                Ts_40.mean(), Ts_50.mean(), Ts_60.mean(), Ts_70.mean(),
#                Ts_80.mean(), Ts_90.mean()])

# print("Standartabweichung Messung 3 {}".format(Ts_std))

# print("**********************************************************\n\n")
# """
# Plots
# """

# p2 = np.polyfit(np.cos(phi)**2, Ts**2, 1)
# x_phi = np.linspace(0, np.radians(90), 100)
# poly2 = np.poly1d(p2)
# Ts_perf = poly2(np.cos(x_phi)**2)
# dTs = np.sqrt((2*Ts*0.4/5)**2)/np.sqrt(3)
# dphi = np.sqrt((-2*np.cos(phi)*np.sin(phi)*np.radians(0.5))**2)




# plt.figure(figsize=(6.4,4))
# #plt.fill_between(d_err**2, Oe_err - dOe_err, Oe_err + dOe_err, color='red', alpha=0.1, label=r'Approximated errorband of $\Theta_m$')
# #plt.plot( np.cos(phi)**2, Ts**2, ".",label=r' $\Theta_c$')
# plt.plot( np.cos(x_phi)**2, Ts_perf, label=r'Approximated linear model')
# #plt.plot(d_err**2, Oe_err, label=r'Calculated $\Theta_c$')
# #plt.plot(d**2, Oe, "or", label=r'Measured $\Theta_m$')
# plt.errorbar(np.cos(phi)**2, Ts**2, xerr=dphi, yerr=dTs, fmt=".", color ="red", capsize=3, label=r'Measured period times $T_m$')
# plt.legend()
# plt.ylabel(r'Time for one period (s)')
# plt.xlabel(r'Alignement of the angle $\cos^2(\varphi)$ ')
# plt.savefig("ellips.pgf")

# # plt.figure()
# # plt.plot(Ts**2, np.cos(phi)**2)
# # plt.show()
