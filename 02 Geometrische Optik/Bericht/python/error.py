# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 14:11:46 2021

@author: fritz
"""
import numpy as np
import matplotlib.pyplot as plt


#Fehlers
def f_obj(a, b):
    return a*b/(a+b)

def del_f_obj(a, b, del_a, del_b):
    dfa = b**2/(a+b)**2
    dfb = a**2/(a+b)**2
    return np.sqrt((dfa*del_a)**2 + (dfb*del_b)**2)



#Messungen:


#Experiment 1


#Objekt- und Bilddistanz
#Linse 1
a11 = 69 #+- 10mm
del_a11 = 0
b11 = 361 #+- 10mm
del_b11 = 10

#calculations
f11 = f_obj(a11, b11)
del_f11 = del_f_obj(a11, b11, del_a11, del_b11)


a12 = 64 #+- 10mm
b12 = 666 #+- 10mm
del_a12 = 10
del_b12 = 0

f12 = f_obj(a12, b12)
del_f12 = del_f_obj(a12, b12, del_a12, del_b12)

#Linse 2
a21=239.5 #+- 10mm
b21=390.5 #+- 10mm

a22=208 #+- 10mm
b22=522 #+- 10mm

#Bessel
#Linse 1
e1 = 67 #+- 10mm
d1 = 240 #+- 10mm

#Linse 2
e2 = 86 #+- 10mm
d2 = 610 #+- 10mm

#---------------------------------------------
#Experiment 2
#Messung 1
a1 = 291 #+- 10mm
b1 = 739 #+- 10mm

#Messung 2
a2 = 332 #+- 10mm
b2 = 568 #+- 10mm


#---------------------------------------------
#Experiment 3
#grob
a_3 = 183 #+- 10mm
b_3 = 857 #+- 10mm
g_3 = 6.8

#fein
a_4 = 176 #+- 10mm
b_4 = 964 #+- 10mm
g_4 = 3.9

#---------------------------------------------
#Experiment 4

lam = 525e-6

#grob
a_5 = 58 #+- 10mm
b_5 = 2842
d_5 = 2.8

#fein
a_6 = 58
b_6 = 2958
d_6 = 2.7
