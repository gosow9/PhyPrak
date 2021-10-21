# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 19:02:43 2021

@author: fritz
"""

#bessel method
def get_f_bessel(e, d):
    return (d**2-e**2)/(4*d)

def lens_eq(a, b):
    return a*b/(a+b)

def d_ber(f, lam, g):
    return f*lam/g

def zers(a, b):
    return 0

def get_div(f_l, f_gem):
    return f_l*f_gem/(f_l - f_gem)

def get_v(a, b):
    return b/a

#experiment 1
a11 = 69
b11 = 361

a12 = 64
b12 = 666

a21=239.5
b21=390.5

a22=208
b22=522



f11 = lens_eq(a11, b11)
f12 = lens_eq(a12, b12)

f_l1 = 0.5*(f11 + f12)

f21 = lens_eq(a21, b21)
f22 = lens_eq(a22, b22)

f_l2 = 0.5*(f21 + f22)

#---------------------------------------------
#experiment 2

e1 = 67
d1 = 240
f_bess_1 = get_f_bessel(e1, d1)


e2 = 86
d2 = 610
f_bess_2 = get_f_bessel(e2, d2)

#---------------------------------------------
#experiment 2
"""
Die Messig wemmer ned, blödi Linse
a1 = 84
b1 = 216

#a1 = 68
#b1 = 632
f_1 = lens_eq(a1, b1)

f_div_1 = get_div(f_l1, f_1)
"""
#messung 1
#a2 = 291
#b2 = 739

#messung 2
a2 = 332
b2 = 568
f_2 = lens_eq(a2, b2)


f_div_2 = get_div(f_l2, f_2)

#---------------------------------------------
#experiment 3

#grob
a_3 = 183
b_3 = 857
g_3 = 6.8
v_3 = get_v(a_3, b_3)


#fein
a_4 = 176
b_4 = 964
g_4 = 3.9
v_4 = get_v(a_4, b_4)

g_grob = g_3/v_3
g_fein = g_4/v_4

#---------------------------------------------
#experiment 4
g = g_grob
lam = 525e-6 #nm
f = f_l2

print(f, lam, g)
d_calc = f*lam/g

a_5 = 58
b_5 = 2842
d_5 = 2.8
v_5 = get_v(a_5, b_5)

a_6 = 58
b_6 = 2958
d_6 = 2.7
v_6 = get_v(a_6, b_6)

d_spalt_grob = d_5/v_5
d_spalt_fein = d_6/v_6





