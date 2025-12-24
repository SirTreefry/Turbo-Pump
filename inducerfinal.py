# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 17:19:01 2025

@author: trefr
"""
import numpy as np
import math
Q = 0.796  # volumetric flow rate [m^3/s]
N = 35000          # rotational speed [RPM]
rho = 71.2          # fluid density [kg/m^3] (LOX example)
g = 32.2           # gravity [m/s^2]

NPSH_a = 39
NPSH_c = NPSH_a / 1.3

g = 32.2
A = .007
H = 300
#geometry
#=====================
omega = 2*np.pi*N/60


L_d_rat = 0.4



#inlet Cm axial velocity
hub_tip_ratio = 0.4


dt = (A/((np.pi/4)*(1 - hub_tip_ratio**2)))**.5

dhub = dt*hub_tip_ratio

Li = dt * L_d_rat

#Ut mean tip speed
Utip = (N*2*dt*np.pi)/60


#Inlet conditions
#==============================
angle_tip = 7*(np.pi/180)#tip half angle
angle_hub = 14*(np.pi/180)#hub half angle
vane_angle_tip = 8


dt0 = dt + Li*math.tan(angle_tip)

dh0 = dhub - Li*math.tan(angle_hub)

Cm0 = Q/(3.12*(np.pi/4)*(dt0**2 - dh0**2))

#inducer_effectd 
d0 = ((dt0**2 + dh0**2)/2)**.5


#inlet peripheral velcotiy
u0 = (N*2*d0*np.pi)/60
#inlet angle beta
B0 = math.atan(Cm0/u0)*(180/np.pi)

uot0 = (N*2*dt0*np.pi)/60

#realtive tip flow angle
B0t = math.atan(Cm0/uot0)*(180/np.pi)
vane0 = dt0/d0 * math.tan(vane_angle_tip*(np.pi/180))*(180/np.pi)
vanle_angle_hub0 = dt0/dh0 * math.tan(vane_angle_tip*(np.pi/180))*(180/np.pi)
#outlet conditions
#==============================
dt1 = dt - Li*math.tan(angle_tip)

dh1 = dhub + Li*math.tan(angle_hub)

Cm1 = Q/(3.12*(np.pi/4)*(dt1**2 - dh1**2))


#inducer_effectd 
d1 = ((dt0**2 + dh0**2)/2)**.5


#inlet peripheral velcotiy
u1 = (N*2*d0*np.pi)/60

cu1 = H*(g/u1)


uot1 = (N*2*dt1*np.pi)/60

#inducer outlet absolute angle
a1 = math.atan(Cm1/cu1)*(180/np.pi)

#inlet angle beta
B1 = math.atan(Cm1/(u1 - cu1))*(180/np.pi)

vane1 = 14.3

vane_angle_tip1 = d1/dt1 * math.tan(vane1*(np.pi/180))*(180/np.pi)

vane_angle_hub1 = d1/dh1 * math.tan(vane1*(np.pi/180))*(180/np.pi)



#Generally geometry
#================================
z = 3


#vane ptich at mean tip diameter
Pi = (np.pi *dt)/z

#chord length at vane tip inlet
Ci = Li / math.sin((vane_angle_tip + vane_angle_tip1)*(np.pi/180)/2)

Solidity = Ci/Pi






#Cavitation Check
#==================================
phi = Cm0/uot0

Ss = N*(Q*60)**.5 *NPSH_c**(-.75)

Ss_correct = Ss/(1 - hub_tip_ratio**2)**.5

phi_check = (3574/Ss_correct)/((1 + (1 + 6*(3574/Ss_correct)**2)**.5)/2)

K = (2*phi**2)/(1 - 2*phi**2)

tau = K + K*phi**2 + phi**2

q_correct = (Q*60)/(1 - hub_tip_ratio**2)


NPSH_r = ((N*q_correct**.5)/Ss_correct)**(4/3)

Ds = d0 * (NPSH_a)**.25 * (q_correct)**(-.5)



print(f"Phi: {phi:.4f}")
print(f"K: {K:.4f}")
print(f"Suction Speed max: {Ss:.2f}")
print(f"Suction Speed Corrected max: {Ss_correct:.2f}")
print(f"NPSHa: {NPSH_a:.4f}")
print(f"NPSHR: {NPSH_r:.4f}")
print(f"Tau: {tau:.4f}")

print("\nInlet Conditions")
print(f"Diameter: {dt0:.4f}")
print(f"Hub diameter : {dh0:.4f}")
print(f"Merirdional inlet velocity : {Cm0:.2f}")
print(f"Tip speed : {uot0:.2f}")
print(f"Vane Tip Angle : {vane_angle_tip:.2f}")

print("\nOUtlet Conditions")
print(f"Diameter: {dt1:.4f}")
print(f"Hub diameter : {dh1:.4f}")
print(f"Merirional Outlet Axial velocity : {Cm1:.2f}")
print(f"Tip speed : {uot1:.2f}")
print(f"Vane Tip Angle : {vane_angle_tip1:.2f}")




