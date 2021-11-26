# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 11:33:43 2021

@author: fritz
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from gauss import error
from matplotlib.ticker import LinearLocator
from scipy.optimize import curve_fit, least_squares
from numpy.linalg import lstsq




print('------------------------------------')
print('Experiment 1')
print('------------------------------------')

V_1 = np.array([200, 300, 400, 500, 600, 700, 800, 900, 980]) #start with 200, +100
t_1 = np.ones_like(V_1)*30#so 1% is reached
N_1 = np.array([0, 0, 1826, 2379, 2535, 2452, 2501, 2402, 2516])
poly = np.polyfit(V_1[3:], N_1[3:], 1)
lin = lambda V: poly[0]*V + poly[1]
V_poly = np.linspace(500, 980)


# plt.figure(figsize=(6.4, 4))
# plt.plot(V_poly, lin(V_poly))
# plt.errorbar(V_1, N_1, None, 10, '-o')
# plt.ylabel('Number of decays, N')
# plt.xlabel('Voltage, V')
# plt.show()

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
N_bg = 50
Neff_2 = N_2-N_bg

d = 101 #mm
r = 18.05 #mm Radius Stahlblende #pm 0.025
theta = (np.tan(r/2/d)) #deg Halbwinkel

eps = np.sin(theta/2)**2

aktiv_2 = Neff_2/t_2/eps #Bq
aktiv_tilde = aktiv_2/2*(100/45+100/90) #Bq
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


print('Halbwinkel: ', np.degrees(theta), 'deg')
print('Akzeptanz: ', eps)
print('Aktivität: ', aktiv_2, 'Bq')
print('Aktivität ohne Folie: ', aktiv_tilde, 'Bq')


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
#----------------------------------------------------------------Cedric-----------------------------------------------------------------------------##

print('\n------------------------------------')
print('Experiment 4')
print('------------------------------------')

# measured values
x_mm = np.array([0, 0.3, 0.5, 1.1, 1.6, 2.,  2.3, 2.5, 2.8, 3., 3.5, 4.])  # Plate thickness in mm
t_3 = 90
N3_tot = np.array([7390, 4606, 3124, 1496, 689, 308, 205, 128, 82, 67, 51, 53])
N3_err = np.sqrt(N3_tot)
N3_sec = N3_tot/t_3

# Change of measurement
x_cm = x_mm/10  # x in cm
n = 1
rho_g = 2.69 #mg cm^-3
rho_mg = 2.69*1000 #mg cm^-3
print("Dichte mal Rmax", rho_mg*x_cm)
N0 = 7390
fun = lambda u, x, y: N0*np.exp(-x*u)+(N_bg) - y

res = least_squares(fun, 1.0, args=( x_mm[:-2], N3_tot[:-2]))
mu = res["x"]
print("fehler vektor:", np.round(np.sqrt(N3_tot)))

N3_log = np.log(N3_tot[:-2])

res1 = np.polyfit(x_mm[:-2], N3_log, 1, full=True)
mu = res1[0]
err = res1[1]

E_max1 = 1.87
E_max2 = (17*rho_g/abs(mu[0])/10)**(1/1.43)
dE_max2 = 100*17**(100/143)*rho_g**(100/143)/(143*abs(mu[0])**(243/143))*err
print("E_max abgelesen = ", E_max1,"MeV" )
print("E_max quelle = ", E_max1 +0.108, "MeV" )
print("E_max berechnet = ", E_max2, "MeV +- ",dE_max2 )
print("E_max Quelle berechnet = ", E_max2+0.11, "MeV+- ",dE_max2+0.2 )
print("E_max Literatur = ", 2.282, "MeV" )
print(17*rho_g*E_max2**(-1.43)/10)


# Plot 1. Methode E-max
N_3 = N3_tot
x_plt = np.linspace(0, 3.2, 100)
plt.figure(figsize=(6.4, 4))
plt.plot(x_plt, N0*np.exp(mu[0]*x_plt), "--",alpha = 0.2, color= "black", label= r'Fitted $N_{x_0} e^{-\mu x}$ , $\mu =$' + "{:.3}".format(abs(mu[0]))+r'$\pm$'+"{:.1}".format(err[0]))
plt.errorbar(x_mm, N3_tot, yerr=N3_err, fmt="o", color ="black", capsize=3, markersize=3,elinewidth=1,fillstyle='none', label=r'Counts per 90 seconds N $\pm \sqrt{N}$')
#plt.plot(x_mm, N_3, "o", label="Counts per 90 seconds N", fillstyle='none', color= "black")


plt.axhline(y=(N_bg), color="black", linestyle=":", label=r'Background noise $N_{BG}$')
plt.yscale('log')
plt.xlabel(r'$x$ [mm]')
plt.ylabel('N/90s')
plt.legend(fontsize = 10)
plt.savefig("fit.png", dpi=600)


E_x = np.array([1.87, 1.57, 1.47, 1.07, 0.83, 0.70, 0.58, 0.5, 0.37, 0.27])
E_max = np.array([0, 0.3, 0.4, 0.8, 1.04, 1.17, 1.29, 1.37, 1.5, 1.6])
N_y = np.array([7340, 4556, 3074, 1446, 639, 258, 155, 78, 32, 17])

E_logx = np.log(E_x)
N_logy = np.log(N_y)

resi = np.polyfit(E_logx, N_logy, 1, full=True)
poly = resi[0]
err1 = resi[1]
print("Squared residuals:", resi[0])
E_plot = np.linspace(E_logx[0], E_logx[-1], 100)

#----------------------------------------------------------------Plot LogLog-----------------------------------------------------------------------------##
plt.figure()
plt.yscale("log")
plt.xscale("log")
plt.errorbar(np.exp(E_logx), np.exp(N_logy),yerr=np.sqrt(N_y),fmt="o", color ="black", capsize=3, markersize=3,elinewidth=1,fillstyle='none',label=r'$N_{eff} \pm \sqrt{N}$' )
plt.plot(np.exp(E_plot), np.exp(poly[1]+poly[0]*E_plot), "--", alpha=0.3, color="grey",  label=("Slope n = {:.2}".format(poly[0]) + r'$\pm$'+" {:.1}".format(err1[0])))
plt.legend()
plt.xlabel(r'$\ln(E_{max}-E(x))$')
plt.ylabel('$\ln(N_{eff})$')
plt.savefig("loglog.png", dpi=600)

#----------------------------------------------------------------Plot Zweiti stiegig-----------------------------------------------------------------------------##
res3 = np.polyfit(E_max, N_y**(1/poly[0]), 1, full=True)
poly3 = res3[0]
err3 = res3[1]
E_max3 = -poly3[1]/poly3[0]


E_plot2 = np.linspace(0, 2.5, 100)
plt.figure()
plt.plot(E_max, N_y**(1/poly[0]), "o", fillstyle='none', color ="black", label=r'$(N_{eff})^{1/n}$'+", n = {:.2}".format(poly[0])+ r'$\pm$'+" {:.1}".format(err1[0]))
plt.plot(E_plot2, poly3[1]+poly3[0]*E_plot2, "--", alpha=0.3, color="grey",label=(r'f(E)='+" {:.3}".format(poly3[1]) +"{:.3}".format(poly3[0])+r'$\cdot E$'))
plt.errorbar(E_max3, 0,xerr= 0.2, fmt="x", color ="black", capsize=5, markersize=6,elinewidth=2, label=(r'f(0)='+" {:.2}".format(E_max3)+r'$\pm$'+"0.3"))
plt.ylim(bottom=0)
plt.legend()
plt.xlabel(r'$E$ [MeV]')
plt.ylabel('$(N_{eff})^{1/n}$')
plt.savefig("emax.png", dpi=600)
print("E_max3: ", E_max3, " +- ",-poly3[1]/poly3[0]**2*err3[0])

#----------------------------------------------------------------Log Diagramm-----------------------------------------------------------------------------##
N3_eff = N3_tot - N_bg
b = (2000000**(-0.4) * 1000000)**(1/0.6)
a = np.log(2000000/b)/1000
#print(b*np.exp(rho_mg*x_cm*a)/1000000)















