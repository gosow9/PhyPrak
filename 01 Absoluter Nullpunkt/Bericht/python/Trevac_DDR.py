# -*- coding: utf-8 -*-
"""
Created on Wed May 26 17:59:06 2021

@author: sebes
"""
import numpy as np
from scipy import stats as stats
import matplotlib.pyplot as plt
from scipy.optimize import leastsq, fsolve
import scipy.optimize as optimization
from scipy.stats import t
import scipy
import pandas as pd
import tikzplotlib

"""
reference values
--> reference sources in report 
"""
#molar masses
molarH = 86.18
molarM = 32.04 
molarE = 46.07

#density 
densityH = 0.659 # g/cm^3 at 25 C
densityM = 0.791 # g/cm^3 at 25 C
densityE =  	0.789 # g/cm^3 at 25 C

#vapour enthalpy
enthalpyM = 37.43 # kj/mol at 25 C from meister script page 25

"""
plotting of the cooling curves
"""

dataH = pd.read_csv('n_hexane.txt',sep="  ",names=['time','temperature'])
dataM = pd.read_table('methanole.dat',sep="  ",names=['time','temperature'])
dataE = pd.read_table('ethanol.txt',sep="  ",names=['time','temperature'])

xH = dataH.loc[2220:3075,'time'].values - 444
yH = dataH.loc[2220:3075,'temperature'].values 

#Integrationsbereiche
timeH1 = dataH.loc[2365:2521,'time'].values - 444
thetaH1 = dataH.loc[2365:2521,'temperature'].values 
timeH2 = xH[413:564]
thetaH2 = yH[413:564]
timeH3 = xH[616:778]
thetaH3 = yH[616:778]

#peakflächen n-hexane
peakH1 = np.trapz(timeH1,thetaH1)
peakH2 = np.trapz(timeH2,thetaH2)
peakH3 = np.trapz(timeH3,thetaH3)
peakH = np.array([peakH1,peakH2,peakH3])
peakHmean = np.mean(peakH)
peakHstd = np.std(peakH)

xE = dataE.loc[900:1555,'time'].values - 180
yE = dataE.loc[900:1555,'temperature'].values

#Integrationsbereiche
timeE1 = xE[25:227]
thetaE1 = yE[25:227]
timeE2 = xE[230:443]
thetaE2 = yE[230:443]
timeE3 = xE[455:645]
thetaE3 = yE[455:645]

#peakflächen ethanol
peakE1 = np.trapz(timeE1,thetaE1)
peakE2 = np.trapz(timeE2,thetaE2)
peakE3 = np.trapz(timeE3,thetaE3)
peakE = np.array([peakE1,peakE2,peakE3])
peakEmean = np.mean(peakE)
peakEstd = np.std(peakE)

xM = dataM.loc[200:1159,'time'].values - 40
yM = dataM.loc[200:1159,'temperature'].values

#Integrationsbereiche
timeM1 = xM[126:332]
thetaM1 = yM[126:332]
timeM2 = xM[403:620]
thetaM2 = yM[403:620]
timeM3 = xM[678:898]
thetaM3 = yM[678:898]

#peakflächen methanol
peakM1 = np.trapz(timeM1,thetaM1)
peakM2 = np.trapz(timeM2,thetaM2)
peakM3 = np.trapz(timeM3,thetaM3)
peakM = np.array([peakM1,peakM2,peakM3])
peakMmean = np.mean(peakM)
peakMstd = np.std(peakM)

"""
Definition for the vapour enthalpy with methanol as reference
"""

def ent(peak,molar,rho):
    
    """
    molar = molar mass from the liquid
    rho = rho from the liquid
    peak = area under peak 

    """
    return (enthalpyM*densityM*molar*peak)/(molarM*rho*peakMmean)

"""
vapour enthalpy calculations and error calculations
"""
enthalpyH = ent(peakHmean,molarH,densityH)
enthalpyE = ent(peakEmean,molarE,densityE)

ts = stats.t(df=2).ppf(0.975)

enthalpystdE = np.sqrt((peakEstd/peakEmean)**2 + (peakMstd/peakMmean)**2)*enthalpyE
verintE = ts*enthalpystdE

enthalpystdH = np.sqrt((peakHstd/peakHmean)**2 + (peakMstd/peakMmean)**2)*enthalpyH
verintH = ts*enthalpystdH

plt.plot(xH,yH)
plt.plot(timeH1,thetaH1,c = 'k')
plt.axvline(x = timeH1[0],linewidth=1, color='k',linestyle='--')
plt.axvline(x = timeH1[156],linewidth=1, color='k',linestyle='--')
plt.fill_between(timeH1, thetaH1[0], thetaH1,color='blue')
plt.plot(timeH2,thetaH2 ,c ='k')
plt.axvline(x = timeH2[0],linewidth=1, color='k',linestyle='--')
plt.axvline(x = timeH2[150],linewidth=1, color='k',linestyle='--')
plt.fill_between(timeH2, thetaH2[0], thetaH2,color='blue')
plt.plot(timeH3,thetaH3 ,c ='k')
plt.axvline(x = timeH3[0],linewidth=1, color='k',linestyle='--')
plt.axvline(x = timeH3[161],linewidth=1, color='k',linestyle='--')
plt.fill_between(timeH3, thetaH3[0], thetaH3,color='blue')
plt.ylabel(r'time / s')
plt.xlabel(r'temperature $\theta$ / C$^{circle}$')
plt.show()
tikzplotlib.save("DDR3.tex")
plt.figure()
plt.plot(xM,yM)
plt.plot(timeM1,thetaM1,c = 'k')
plt.axvline(x = timeM1[0],linewidth=1, color='k',linestyle='--')
plt.axvline(x = timeM1[205],linewidth=1, color='k',linestyle='--')
plt.fill_between(timeM1, thetaM1[0], thetaM1,color='blue')
plt.plot(timeM2,thetaM2,c = 'k')
plt.axvline(x = timeM2[0],linewidth=1, color='k',linestyle='--')
plt.axvline(x = timeM2[216],linewidth=1, color='k',linestyle='--')
plt.fill_between(timeM2, thetaM2[0], thetaM2,color='blue')
plt.plot(timeM3,thetaM3,c = 'k')
plt.axvline(x = timeM3[0],linewidth=1, color='k',linestyle='--')
plt.axvline(x = timeM3[219],linewidth=1, color='k',linestyle='--')
plt.fill_between(timeM3, thetaM3[0], thetaM3,color='blue')
plt.ylabel(r'time / s')
plt.xlabel(r'temperature $\theta$ / C$^{circle}$')
plt.show()
tikzplotlib.save("DDR4.tex")
plt.figure()
plt.plot(xE,yE)
plt.plot(timeE1,thetaE1,c = 'k')
plt.axvline(x = timeE1[0],linewidth=1, color='k',linestyle='--')
plt.axvline(x = timeE1[201],linewidth=1, color='k',linestyle='--')
plt.fill_between(timeE1, thetaE1[0], thetaE1,color='blue')
plt.plot(timeE2,thetaE2,c = 'k')
plt.axvline(x = timeE2[0],linewidth=1, color='k',linestyle='--')
plt.axvline(x = timeE2[212],linewidth=1, color='k',linestyle='--')
plt.fill_between(timeE2, thetaE2[0], thetaE2,color='blue')
plt.plot(timeE3,thetaE3,c = 'k')
plt.axvline(x = timeE3[0],linewidth=1, color='k',linestyle='--')
plt.axvline(x = timeE3[189],linewidth=1, color='k',linestyle='--')
plt.fill_between(timeE3, thetaE3[0], thetaE3,color='blue')
plt.ylabel(r'time / s')
plt.xlabel(r'temperature $\theta$ / C$^{circle}$')
plt.show()
tikzplotlib.save("DDR5.tex")