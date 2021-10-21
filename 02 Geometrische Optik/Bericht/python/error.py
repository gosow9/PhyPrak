# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 14:11:46 2021

@author: fritz
"""
import numpy as np


#Fehlers


#Messungen:


#Experiment 1
def f_obj(a, b):
    return a*b/(a+b)

def del_f_obj(a, b, del_a, del_b):
    dfa = b**2/(a+b)**2
    dfb = a**2/(a+b)**2
    return np.sqrt((dfa*del_a)**2 + (dfb*del_b)**2)

#Objekt- und Bilddistanz
#Linse 1
a11 = 69 #+- 5mm
del_a11 = 5
b11 = 361
del_b11 = 0

#calculations
f11 = f_obj(a11, b11)
del_f11 = del_f_obj(a11, b11, del_a11, del_b11)


a12 = 64 #+- 5mm
b12 = 666
del_a12 = 5
del_b12 = 0

f12 = f_obj(a12, b12)
del_f12 = del_f_obj(a12, b12, del_a12, del_b12)

f1_l = 0.5*(f11+f12)
del_f1_l = 0.5*np.sqrt((del_f11**2+del_f12**2))

#Linse 2
a21=239.5 #+- 10mm
b21=390.5

del_a21 = 10
del_b21 = 0

f21 = f_obj(a21, b21)
del_f21 = del_f_obj(a21, b21, del_a21, del_b21)

a22=208 #+- 10mm
b22=522
del_a22 = 10
del_b22 = 0

f22 = f_obj(a22, b22)
del_f22 = del_f_obj(a22, b22, del_a22, del_b22)

f2_l = 0.5*(f21+f22)
del_f2_l = 0.5*np.sqrt((del_f21**2+del_f22**2))

#Bessel

def f_bes(e, d):
    return (d**2 - e**2)/(4*d)

def del_f_bes(e, d, del_e):
    dfe = -0.5*e/d
    return abs(dfe*del_e)


del_e = 20
#Linse 1
e1 = 67 #+- 20mm
d1 = 240 #+- 0.5mm, vernachlässigbar
f1_b = f_bes(e1, d1)
del_f1_b = del_f_bes(e1, d1, del_e)

#Linse 2
e2 = 86 #+- 20mm
d2 = 610 #+- 0.5mm
f2_b = f_bes(e2, d2)
del_f2_b = del_f_bes(e2, d2, del_e)

#---------------------------------------------
#Experiment 2
def f_div(f_conv, f_gem):
    return f_conv*f_gem/(f_conv - f_gem)

def del_f_div(f_conv, f_gem, del_f_conv, del_f_gem):
    dfc = -f_gem**2/(f_conv - f_gem)**2
    dfg = f_conv**2/(f_conv - f_gem)**2
    return np.sqrt((dfc*del_f_conv)**2 + (dfg*del_f_gem)**2)

#Messung 1
a1 = 291 #+- 10mm
b1 = 739 
f1 = f_obj(a1, b1)
del_f1 = del_f_obj(a1, b1, 10, 0)
f_div1 = f_div(f2_b, f1)
del_f_div1 = del_f_div(f2_b, f1, del_f2_b, del_f1)


#Messung 2
a2 = 332 #+- 10mm
b2 = 568 
f2 = f_obj(a2, b2)
del_f2 = del_f_obj(a2, b2, 10, 0)
f_div2 = f_div(f2_b, f2)
del_f_div2 = del_f_div(f2_b, f2, del_f2_b, del_f2)

#---------------------------------------------
#Experiment 3
def del_v(a, b, del_a=10):
    dva = -b/a**2
    return abs(dva*del_a)

def del_g(g_t, v, del_v, del_g_t = 0.1):
    dgg = 1/v
    dgv = -g_t/v**2
    return np.sqrt((dgg*del_g_t)**2 + (dgv*del_v)**2)
    
#grob
a3 = 183 #+- 10mm
b3 = 857 
g3_t = 6.8 #+- 0.1 mm
v3 = b3/a3
del_v3 = del_v(a3, b3)
g3 = g3_t/v3
del_g3 = del_g(g3, v3, del_v3)

#fein
a4 = 176 #+- 10mm
b4 = 964
g4_t = 3.9 #+- 0.1 mm
v4 = b4/a4
del_v4 = del_v(a4, b4)
g4 = g4_t/v4
del_g4 = del_g(g4, v4, del_v4)
#---------------------------------------------
#Experiment 4

lam = 525e-6 #+- 0.5e-6
del_lam = 0.5e-6


#grob
a5 = 58 #+- 10mm
b5 = 2842
v5 = b5/a5
del_v5 = del_v(a5, b5)
d5_t = 2.8
d5 = d5_t/v5
del_d5 = del_g(d5_t, v5, del_v5)


#fein
a6 = 58
b6 = 2958
v6 = b6/a6
del_v6 = del_v(a6, b6)
d6_t = 2.7
d6 = d6_t/v6
del_d6 = del_g(d6_t, v6, del_v6)

#abbe
def del_abbe(f, lam, g, del_f, del_lam, del_g):
    ddf = lam/g
    ddl = f/g
    ddg = -f*lam/g**2
    return np.sqrt((ddf*del_f)**2 + (ddl*del_lam)**2 + (ddg*del_g)**2)

d_calc = f2_b*lam/g3
del_d_calc = del_abbe(f2_b, lam, g3, del_f2_b, del_lam, del_g3)

