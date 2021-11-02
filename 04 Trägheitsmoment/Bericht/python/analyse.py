# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 12:24:54 2021

@author: Cedric
"""

import numpy as np
import matplotlib.pyplot as plt


def moment_disc(Ms, Rs):
    return 0.5*Ms*Rs**2


def torsion_disc(Os, T1, T2):
    return 4*np.pi**2*Os/(T2**2 - T1**2)


def moment_quad(M, b, c):
    return M*(b**2 + c**2)/12


def steiner(Os, d, M):
    return Os + d**2*M


def moment_ta(Ta, D):
    return D*Ta**2/(8*np.pi**2)


"""
Experiment 1 bestimmung von D
"""
print("**********************************************************")
print("Experiment 1")
print("----------------------------------------------------------\n")


T1 = np.array([32.59, 32.34, 32.37, 32.47, 32.53])/5  # Measurements in second
T2 = np.array([25.63, 25.59, 25.50, 25.53, 25.47])/3  # Measurements in second
print("Mean T1 = {}, std of T1 = {}".format(T1.mean(), T1.std()))
print("Mean T2 = {}, std of T2 = {}".format(T2.mean(), T2.std()))
Ms = 3.431             # Disc weight in kg +- 0.5 gramm
Rs = 0.249/2       # Disc radius in m
Os = moment_disc(Ms, Rs)
D = torsion_disc(Os, T1.mean(), T2.mean())
print(r'Phi disc = {}'.format(Os))
print(r'D disc = {}'.format(D))
print("**********************************************************\n\n")


"""
Experiment 2 Satz von Steiner
"""
print("Experiment 2")
print("----------------------------------------------------------\n")


Ta_5 = np.array([20.34, 20.38, 20.40])/3
Ta_10 = np.array([22.50, 22.38, 22.50])/3
Ta_15 = np.array([25.72, 25.69, 25.65])/3
Ta_20 = np.array([29.71, 29.75, 29.69])/3
Ta_25 = np.array([34.15, 33.97, 34.22])/3
Ta = np.array([T1.mean(), Ta_5.mean(), Ta_10.mean(),
               Ta_15.mean(), Ta_20.mean(), Ta_25.mean()])
Mq = np.array([0.611, 0.611, 0.611, 0.605, 0.605, 0.605])  # Quad weight in kg
c = 0.6
d = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25])
b = 0.01174

# ----------------------------------------------------------------------------

print("Mean 1 bar = {}kg, std of bars = {}, Two bars = {}kg".format(
    Mq.mean(), Mq.std(), 2*Mq.mean()))
Oe = moment_ta(Ta, D)             # Experimental moment of inertia
Obar = moment_quad(Mq.mean(), b, c)  # Calculated momentof inertia
Oc = steiner(Obar, d, Mq.mean())
print("Theta experimental = {}".format(Oe))
print("Theta calculated = {}".format(Oc))
print("**********************************************************\n\n")

"""
Experiment 3 Ellipse
"""
print("Experiment 3")
print("----------------------------------------------------------\n")
Ts_0 = np.array([17.28, 17.37, 17.34])/5
Ts_10 = np.array([17.65, 17.63, 17.66])/5
Ts_20 = np.array([18.25, 18.28, 18.31])/5
Ts_30 = np.array([19.34, 19.31, 19.25])/5
Ts_40 = np.array([20.44, 20.40, 20.47])/5
Ts_50 = np.array([21.59, 21.59, 21.53])/5
Ts_60 = np.array([22.57, 22.46, 22.50])/5
Ts_70 = np.array([23.31, 23.31, 23.35])/5
Ts_80 = np.array([23.97, 23.85, 23.84])/5
Ts_90 = np.array([23.97, 24.03, 24.06])/5
phi = np.radians(np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90]))  # 0.2 Grad error
print(np.cos(phi)**2)

Ts_std = np.array([Ts_0.std(), Ts_10.std(), Ts_20.std(), Ts_30.std(),
               Ts_40.std(), Ts_50.std(), Ts_60.std(), Ts_70.std(),
               Ts_80.std(), Ts_90.std()])
Ts = np.array([Ts_0.mean(), Ts_10.mean(), Ts_20.mean(), Ts_30.mean(),
               Ts_40.mean(), Ts_50.mean(), Ts_60.mean(), Ts_70.mean(),
               Ts_80.mean(), Ts_90.mean()])

print("Standartabweichung Messung 3 {}".format(Ts_std))

print("**********************************************************\n\n")
"""
Plots
"""
plt.figure()
plt.plot(Oe, d**2)
plt.show()

plt.figure()
plt.plot(Ts**2, np.cos(phi)**2)
plt.show()