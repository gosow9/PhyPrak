import matplotlib.pyplot as plt
import numpy as np




u1 = 132.21
u2 = -3.72
p1 = 95711.7612 
p2 = 10 #gegeben, +-10

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
    t0_2 = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
    return t0_1, t0_2

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
    c = (p2-p1)/(u2-u1) #Steigung Druck/Spannung
    #print("d", c)
    p0 = (p2-c*u2) #Offset
    #print("p0 :", p0)
    return p0 + c*U

ue = 96.18
uk = 131.54
un = 24.63

pe = calc_p(ue)
pk = calc_p(uk)
tk = 98.412
tl = 23
te = 0
pn = calc_p(un)

t0 = calc_t0(pe, pk, tk, tl)[0]
t0_ex = -273.15 #exact

print("t0 :", t0)
print("nitrogen: ", calc_tn(pe, pn, te, t0, tl))
