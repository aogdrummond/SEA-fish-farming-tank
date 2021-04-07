# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 23:58:56 2021

@author: aogdrummond
"""

import numpy as np
import math
#%% System's Constants

rho0 = 1.2922        #air density [kg/m³]
c0 = 343.0           #air sound speed [m/s]
ref = 20e-6          #pressure reference [Pa]
rho1 = 1000.0        #water density [kg/m³]
c1 = 1480.0;         #water sound speed [m/s]

#%%Subsystems dimensions

#Rooms dimensions
r1 = 5.65
r2 = 4.1
r3 = 3.2

#Tanks dimensions
t1 = 3.962
t2 = 3.048
t3 = 2.59

Vr = r1*r2*r3        #room's volume [m³]
Vt = t1*t2*t3        #tank's volume [m³]


Sroom =  2*(r1*r2) + 2*(r1*r3) + 2*(r2*r3) - (t1*t3)  #Sum of the areas
Lroom =  4*r1+4*r2+4*r3               # Sum of the sides [m]
Stank =  2*(t1*t2)+2*(t2*t3) + 2*(t1*t3)
Ltank =  4*t1+4*t2+4*t3
Atank = t1*t3      #Tanks Surface area

#%% Reverberation Time

Paper_T60_r = [8.61, 7.54, 7.76, 8.69, 10.43, 10.65, 8.92, 8.59, 7.98, 6.74, 6.07, 5.34, 4.96, 4.41, 3.73, 3.05, 2.68, 2.31];
Paper_T60_t = [0, 0, 0, 0, 0, 0, 0, 0, .35, .42, .43, .33, .32, .25, .26, .24, .25, .25];

#%% Vacuum Sound Power

We = [1.067701213366e-7, 3.258485333e-6, 2.396535e-7, 4.50449e-7, 1.349024e-6, 2.190076e-6, 3.531653e-6, 5.11264e-6, 0.00024915652, 3.21072299e-5, 7.0694245e-6, 7.5555327e-6, 1.1512413e-5, 1.489275e-5, 1.7656767e-5, 8.76629240e-6, 6.7690471e-6, 6.984186e-6];

  
#%% Getting pump's sound power

W0 = 10e-12                   #Reference Power
Power = sum(We)               # Total Power (summing each band)
Original_SWL = 10*np.log10(abs((Power)/W0))

We = np.multiply(We,42)                    #Power Compensation for the power to be 92 dB
Adapted_Power = sum(We)
Adapted_SWL = 10*np.log10(abs((Adapted_Power)/W0))

#%% Schroeder's Frequency

Fsr = 2000*np.sqrt(np.mean(Paper_T60_r)/Vr);    #Rooms Schroeder's frequency
Fst = 2000*np.sqrt(np.mean(Paper_T60_t)/Vt);    #Tank's Schroeder's frequency

#%% Array of Frequencies

freq_bands = [100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000]         #Frequencies (one-third octave bands)


for i in range(len(freq_bands)):
    if freq_bands[i] >= Fsr:
        break

f = freq_bands[(i):]       #Use only frequencies above Schroeder's frequency
w = 2*np.multiply((math.pi),f)              #From linear to angular freq.

T60_r = Paper_T60_r[(i):]
T60_t = Paper_T60_t[(i):]


#%% Energies Calculation

Energy_Mat = []
for i in range(len(f)):
    
    #Losses Factors
    
    etar= (2.2/(f[i]*T60_r[i]))
    m_t = (rho1*Vt)/(t1*t3)        
    a_t = (m_t*t1*w[i])/(2*rho1*c1) 
    tau_rt = np.log(1+(a_t**2))/(a_t**2)
    etart= ((c1*Atank)/(4*w[i]*Vr))*tau_rt
    etat= (2.2/(f[i]*T60_t[i]))
    m_r = (rho0*Vr)/(t1*t3)
    a_r = (m_r*t3*w[i])/(2*rho0*c0)
    tau_tr = np.log(1+(a_r**2))/(a_r**2)
    etatr= ((c0*Atank)/(4*w[i]*Vt))*tau_tr
    
    
    ETA = [[-etart, (etat + etatr)],[(etar + etart),-etatr]]   #Coefficients Matrix
    W = [0,We[i]]                                              #Injected Power Matrix
    
    Energy_Mat.append(np.multiply(W,np.asarray(np.multiply(w[i],ETA))))
    
    
#%% Spatial Average Sound Pressure

Pref1 = 2e-5                                  #Air Sound pressure reference
Pref2 = 1e-6                                  #Underwater sound pressure reference
Pr = np.sqrt(Energy_Mat[1]*rho0*(c0**2)/Vr)   #Average Sound Pressure inside the room
Pt = np.sqrt(Energy_Mat[2]*rho1*(c1**2)/Vt)   #Average Sound Pressure inside the tank

#To dB
Pr_dB = 20*np.log10(abs(Pr/Pref1))
Pt_dB = 20*np.log10(abs(Pt/Pref2))

#%% Plots

import matplotlib.pyplot as plt

mean_SPL = np.mean(Pr_dB)
mean_SPL = mean_SPL * np.ones((len(f),1))
total_SPL= 20*np.log10(abs(sum(Pt)/Pref1))

plt.figure()
plt.semilogx(f,Pr_dB)
# txt = ['Spatial Average SPL inside the tank (TOTAL = ' num2str(total_SPL,'%0.0f') 'dB)']
plt.title()
"""
Formatação do gráfico
"""

mean_SPL = np.mean(Pt_dB)
mean_SPL = mean_SPL * np.ones((len(f),1))
total_SPL= 20*np.log10(abs(sum(Pt)/Pref2))

plt.figure()
plt.semilogx(f,Pr_dB)
# txt = ['Spatial Average SPL inside the tank (TOTAL = ' num2str(total_SPL,'%0.0f') 'dB)']
plt.title()
"""
Formatação do gráfico
"""

#%% Sound Source Power

W0 = 10e-12
We_NWS = 10*np.log10(abs((np.asarray(We))/W0))

plt.figure()
plt.semilogx(freq_bands,We_NWS)
plt.title("Pump's Sound Power")

"""
Formatação do gráfico
"""

#%% Room's Reverberation Time

plt.figure()
plt.semilogx(freq_bands, Paper_T60_r)
plt.title("Room's Reverberation Time")
"""
Formatação do gráfico
"""

plt.figure()
plt.semilogx(freq_bands, Paper_T60_t)
plt.title("Tank's Reverberation Time")
"""
Formatação do gráfico
"""

