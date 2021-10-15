# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 15:17:37 2021

@author: fritz
"""

import matplotlib.pyplot as plt
import numpy as np




ul = 131.21
ut = -3.72
pl = 95711.7612 
pt = 10 #gegeben, +-10

def calc_t0(p_e, p_k, t_k, t_l, eps = 0.001, gamma = 10e-5):
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
    #print("a =", a)
    b = eps*(p_k-p_e)*t_k + (1+gamma*t_k)*p_k*t_l - p_e*(t_l + t_k)
    #print("b = ", b)
    c = p_e*t_l*t_k
    #print("c = ", c)
    t0_1 = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    #t0_2 = (-b - np.sqrt(b**2 - 4*a*c))/(2*a) nicht physikalisch
    return t0_1

def calc_tn(pe, pn, te, t0, tl, eps = 0.001, gamma = 10e-5):
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

def calc_p(U):
    """
    Berechnet Druck aus Spannung

    Parameters
    ----------
    U : Spannung in mV

    Returns
    -------
    Druck in Pascal

    """
    c = (pt-pl)/(ut-ul) #Steigung Druck/Spannung
    #print("d", c)
    p0 = (pt-c*ut) #Offset
    #print("p0 :", p0)
    return p0 + c*U

ue = 96.18
uk = 131.28
un = 24.63

pe = calc_p(ue)
pk = calc_p(uk)
tk = 98.412
tl = 23
te = 0
pn = calc_p(un)
c = (pt-pl)/(ut-ul)
p0 = (pt-c*ut)

t0 = calc_t0(pe, pk, tk, tl)
t0_ex = -273.15 #exact



"""
Fehlerberechnung mit Gauss'scher Fehlerfortpflanzung
"""

def find_delta_u(u): #absoluter Messfehler + Fehler in Steigung
    return 0.155+u*0.001




delta_pt = 10
delta_pl = 0.
delta_ut = 0.155
delta_ul = find_delta_u(ul)
delta_tk = 0.
delta_tl = 0.05
delta_te = 0.



A = 1/(ut-ul)*delta_pt
B = -1/(ut-ul)*delta_pl
C = (pl-pt)/(ut-ul)**2*delta_ut
D = (pt-pl)/(ut-ul)**2*delta_ul

delta_c = np.sqrt(A**2 + B**2 + C**2 + D**2)

E = delta_pt
F = -ut*delta_c
G = -c*delta_ut

delta_p0 = np.sqrt(E**2+F**2+G**2)

def find_delta_p(u):
    H = delta_p0
    I = u*delta_c
    J = c*find_delta_u(u)
    return np.sqrt(H**2+I**2+J**2)


delta_pe = find_delta_p(ue)
delta_pk = find_delta_p(uk)


h = 1e-5
dt0_pe = (calc_t0(pe+h, pk, tk, tl) - t0)/h

dt0_pk = (calc_t0(pe, pk+h, tk, tl) - t0)/h

dt0_tk = (calc_t0(pe, pk, tk+h, tl) - t0)/h

dt0_tl = (calc_t0(pe, pk, tk, tl+h) - t0)/h

delta_t0 = np.sqrt((dt0_pe*delta_pe)**2 + (dt0_pk*delta_pk)**2 + (dt0_tk*delta_tk)**2 + (dt0_tl*delta_tl)**2)

tn = calc_tn(pe, pn, te, t0, tl)

delta_pn = find_delta_p(un)

dtn_pe = (calc_tn(pe+h, pn, te, t0, tl)-tn)/h
dtn_pn = (calc_tn(pe, pn+h, te, t0, tl)-tn)/h
dtn_te = (calc_tn(pe, pn, te+h, t0, tl)-tn)/h
dtn_t0 = (calc_tn(pe, pn, te, t0+h, tl)-tn)/h
dtn_tl = (calc_tn(pe, pn, te, t0, tl+h)-tn)/h

delta_tn = np.sqrt((dtn_pe*delta_pe)**2 + (dtn_pn*delta_pn)**2 + (dtn_te*delta_te)**2 + (dtn_t0*delta_t0)**2 + (dtn_tl*delta_tl)**2)


print("t0 :", t0, "+-", delta_t0, "°C")
print("nitrogen: ", calc_tn(pe, pn, te, t0, tl), "+-", delta_tn, "°C")



def plot():

    
    P = np.array([pt, p0, pl])
    U = np.array([ut, 0, ul])
    plt.figure(figsize=(12,9))
    plt.grid()
    plt.xlabel('U [mV]')
    plt.ylabel('p [Pa]')
    plt.title("Calibration")
    plt.plot(U, P, '-')
    plt.errorbar(ul, pl, delta_pl, delta_ul, 'o', label = 'Labor')
    plt.errorbar(ut, pt, delta_pt, delta_ut, 'o', label = 'Vacuum')
    plt.errorbar(0, p0, delta_p0, 0, 'o', label = 'Calculated Offset')
    plt.legend()
    plt.savefig(fname = "Calibration")
    plt.show()
    
    T = np.array([tk, te, t0])
    P_exp = np.array([pk, pe, 0])

    
    
    plt.figure(figsize=(12,9))
    plt.grid()
    plt.xlabel('T [$\degree$ C]')
    plt.ylabel('p [Pa]')
    plt.title("Calculation Absolute Zero")
    plt.plot(T, P_exp, '-')
    plt.errorbar(te, pe, delta_pe, delta_te, 'o', label = 'Ice')
    plt.errorbar(tk, pk, delta_pk, delta_tk, 'o', label = 'Vapor')
    plt.errorbar(t0, 0, 0, delta_t0, 'x', label = 'Absolute Zero')
    plt.errorbar(t0_ex, 0, 0, 0.005, 'x', label = 'Exact value of absolute zero')
    plt.legend()
    plt.savefig(fname = "Calculation Absolute Zero")
    plt.show()
    
    plt.figure(figsize=(12, 9))
    plt.errorbar(t0, 0, 0, delta_t0, fmt='o', label = 'Absolute Zero')
    plt.errorbar(t0_ex, 0, 0, 0.005, 'o', label = 'Exact value of absolute zero')
    plt.xlabel('T [$\degree$ C]')
    plt.ylabel('p [Pa]')
    plt.title("Absolute Zero")
    plt.legend()
    plt.savefig(fname = "Absolute Zero")
    plt.show()
    

plot()







