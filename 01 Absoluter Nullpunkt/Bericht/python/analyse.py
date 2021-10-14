# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 16:29:18 2021

@author: Cedric
"""
import numpy as np
from scipy import stats as stats
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
from scipy.optimize import leastsq, fsolve
import scipy.optimize as optimization
from scipy.stats import t
import scipy
import pandas as pd
import tikzplotlib


# Temperaturen
ts = np.linspace(-280, 110, 1000)
t_luft = 23.1
t_eiswasser = 0
t_Wkochend = 98.412

# Gemessener druck in Pascal
p_luft = 95711.7612
p_vakuum = 10

# Gemessener druck in milivolt
us = np.linspace(-10, 150, 1000)
u_luft = 132.21
u_vakuum = -3.15
u_Wkoch = 131.28


# Sensorcharakteristik

c1 = (p_luft-p_vakuum)/(u_luft-u_vakuum)  # Steigung Druck/Spannung
p1 = (p_luft-c1*u_luft)  # Offset


def calc_t0(p_e, p_k, t_k, t_l, eps=0.001, gamma=10e-5):
    """
    Berechnung absoluter Nullpunkt

    Parameters
    ----------
    p_e : Druck in Eis
    p_k : Druck in Dampf
    t_k : Temperatur Dampf
    t_l : Temperatur Labor

    Returns
    -------
    Beide Lösungen, nur eine physikalisch

    """

    a = (1+eps)*p_e - (1 + eps + gamma*t_k)*p_k
    b = eps*(p_k-p_e)*t_k + (1+gamma*t_k)*p_k*t_l - p_e*(t_l + t_k)
    c = p_e*t_l*t_k
    t0_1 = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    t0_2 = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
    return t0_1, t0_2


def calc_tn(pe, pn, te, t0, tl, eps=0.001, gamma=10e-5):
    """


    Parameters
    ----------
    pe : Druck in Eis
    pn : Druck in Stickstoff
    te : Temperatur Eis
    t0 : Nullpunkt
    tl : Temperatur Luft

    Returns
    -------
    Temperatur Stickstoff

    """
    A = pe/(te-t0) + (eps*(pe-pn))/(tl-t0)
    return (A*t0 + pn)/(A - gamma*pn)


def calc_t00(pe, pk, tk):
    """
    Approximative Berechnung abs. Nullpunkt

    """
    return -pe/(pk-pe)*tk


def p_from_u(U, c, p0):
    """
    Berechnet Druck aus Spannung

    Parameters
    ----------
    U : Spannung in mV

    Returns
    -------
    Druck in Pascal

    """
    return p0 + c*U

error = 0.075
plt.figure()
yti1 = [-20000.,0.,4000,20000.,40000.,60000.,80000.,100000.]
yti2 = [-20000.,0.,2237.1,20000.,40000.,60000.,80000.,100000.]
plt.yticks(yti1,yti2)
#plt.fill_between(us, p_from_u(us-error, c1, p1), p_from_u(us+error, c1, p1))
plt.plot(u_luft, p_luft, 'o', label= r'Ambient pressure')
plt.plot(u_vakuum, p_vakuum, 'o', label= r'Vacuum')
plt.plot(us, p_from_u(us, c1, p1),linestyle='--', linewidth=1, label= r'Calculated slope $c$')
plt.plot(0, p1, 'o', label= r'Calculated offset $p0$')
print(p1)
#plt.vlines(x=0, ymin=-15000, ymax=p1, linewidth=1, color='k', linestyle='--')
#plt.hlines(y=p1, xmin=-20, xmax=0, linewidth=1, color='k', linestyle='--')
plt.xlim(left=-10, right = 150)
plt.ylim(bottom = -1000, top = 100000)
plt.ylabel(r'Pressure in pa')
plt.xlabel(r'Pressure in mV')
plt.legend()


#plt.xticks(list(plt.xticks()[0]) + extraticks)
#plt.grid(True)
plt.savefig("cali.pgf")
#tikzplotlib.clean_figure()
#tikzplotlib.save("cali.tex")
