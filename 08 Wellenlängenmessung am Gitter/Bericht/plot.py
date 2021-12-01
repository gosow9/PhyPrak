# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 15:03:05 2021

@author: Cedric
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from io import StringIO
import matplotlib.transforms as mt


# Webpage url
url = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=hg&limits_type=0&low_w=&upp_w=&unit=1&submit=Retrieve+Data&de=0&format=0&line_out=0&en_unit=0&output=0&bibrefs=1&page_size=15&show_obs_wl=1&show_calc_wl=1&unc_out=1&order_out=0&max_low_enrg=&show_av=2&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&intens_out=on&max_str=&allowed_out=1&forbid_out=1&min_accur=&min_intens=&conf_out=on&term_out=on&enrg_out=on&J_out=on'

# Extract tables
dfs = pd.read_html(url)

# Get first table
df = dfs[3]
df1 = df[df.Ion != 'Hg II' ]
df1 = df1[df1.Ion != 'Hg III' ]
df3 = df1.iloc[:,[1,5]]
df3.rename(columns={ df3.columns[1]: 'int' }, inplace = True)
df3.rename(columns={ df3.columns[0]: 'lam' }, inplace = True)
df3 = df3.dropna()
df3['int'] = df3['int'].map(lambda x: x.lstrip('+-').rstrip('ahHuErcwp,s*'))

# Extract columns
df2 = df.iloc[:,[1,5]]
df2.rename(columns={ df2.columns[1]: 'int' }, inplace = True)
df2.rename(columns={ df2.columns[0]: 'lam' }, inplace = True)
df2 = df2.dropna()
df2['int'] = df2['int'].map(lambda x: x.lstrip('+-').rstrip('ahHuErcwp,s*'))


df_clean = df3[df3.int.str.isnumeric()]
wavelengths = df_clean['lam'].str.replace(' ', '').to_numpy().astype(np.float)[120:250]
spectrum = df_clean['int'].str.replace(' ', '').to_numpy().astype(np.float)[120:250]
spectrum = np.clip(spectrum,0,1000)/1000



def wavelength_to_rgb(wavelength, gamma=0.8):
    ''' taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    Additionally alpha value set to 0.5 outside range
    '''
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 750:
        A = 1.
    else:
        A=0.5
    if wavelength < 380:
        wavelength = 380.
    if wavelength >750:
        wavelength = 750.
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R,G,B,A)

clim=(350,780)
norm = plt.Normalize(*clim)
wl = np.arange(clim[0],clim[1]+1,2)
colorlist = list(zip(norm(wl),[wavelength_to_rgb(w) for w in wl]))
spectralmap = matplotlib.colors.LinearSegmentedColormap.from_list("spectrum", colorlist)

mes1 = np.array([[410.61021899, 442.54749041, 496.49570775, 500.73229807, 551.99395783,
 581.93041616, 584.01148342, 628.39216335, 694.81019667],[410.61021899, 442.54749041, 496.49570775, 500.73229807, 551.99395783,
 581.93041616, 584.01148342, 628.39216335, 694.81019667]])

mes2 = np.array([404.31903822, 435.1152119,  490.7591436,  495.26416482, 545.16959453,
 575.89248023, 578.14537396, 622.56162682, 689.99883605])
ymes = np.ones(9)
ymes2 = np.zeros(9)


fig, axs = plt.subplots(1, 1, figsize=(12,3), tight_layout=True)
# wavelengths = np.linspace(200, 1000, 1000)
# spectrum = (5 + np.sin(wavelengths*0.1)**2) * np.exp(-0.00002*(wavelengths-600)**2)
plt.plot(wavelengths, spectrum, color='black')
plt.plot(mes1,[ymes,ymes2],"--", color='black')

y = np.linspace(0, 0.999, 1000)
X,Y = np.meshgrid(wavelengths, y)

extent=(np.min(wavelengths), np.max(wavelengths), np.min(y), np.max(y))

plt.imshow(X, clim=clim,  extent=extent, cmap=spectralmap, aspect='auto')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity')

plt.fill_between(wavelengths, spectrum, 1.1, color='w')
plt.savefig('WavelengthColors.png')

plt.show()