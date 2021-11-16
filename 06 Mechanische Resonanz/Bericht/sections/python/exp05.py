# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 12:24:54 2021

@author: Cedric
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq, fsolve, curve_fit
import matplotlib
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


alpha1, pcov1 = curve_fit(lambda t, a: 110*np.exp(-a*t),  PT1,  A1, 0.0)
PT_1 = np.linspace(0, PT1[-1],100)
expon1 = 110*np.exp(-alpha1*PT1)
perr = np.std(A1-expon1)
print(perr)

plt.figure()
arr_err = np.ones(A1.shape)*perr
plt.errorbar(PT1, A1, yerr=arr_err, fmt=".", color ="red", capsize=3, label=r'Measured Times')
plt.plot(PT1, expon1)

alpha2 = curve_fit(lambda t, a: 110*np.exp(-a*t),  PT2,  A2, 0.0)
expon2 = 110*np.exp(-alpha2[0]*PT2)

alpha3 = curve_fit(lambda t, a: 110*np.exp(-a*t),  PT3,  A3, 0.0)
expon3 = 110*np.exp(-alpha3[0]*PT3)
print("alpha1=", alpha1)
print("alpha2=", alpha2)
print("alpha3=", alpha3)

plt.figure()
plt.yscale("log")
plt.plot(PT1, A1/A0, ".", label="Dämpfung 1")
plt.plot(PT1, expon1/A0, label="Model 1 ")
plt.plot(PT2, A2/A0, ".", label="Dämpfung 2")
plt.plot(PT2, expon2/A0, label="Model 2")
plt.plot(PT3, A3/A0, ".", label="Dämpfung 3")
plt.plot(PT3, expon3/A0, label="Model 13")










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
