# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 20:53:49 2026

@author: trefr
"""
import numpy as np
import math
import matplotlib.pyplot as plt

#1D Pump Sizing for impellers
#Based on Pumpas method
g = 32.2

#Inputs dimensions
#===========================================

#number of blades
z = 5

#normal blade thickeness
thk = .15

m = 11.02

mL = .1

rho1 = 62.4

N = 3600

Q = m/rho1

Pt1 = 100

Pv = 50.14

K = 3*10**(-12)

#Entrance or stage 1 values
#=======================================

#blade span hub to tip
B1 = 1.77

#blade angle
beta1 = 20

alpha = 25

Rhub1 = 1.1811

Rtip1 = 2.953

lamba1 = 1 - .01

S = .1

#Entrance or stage 2 values
#=======================================
Rhub2 = 5.9

Rtip2 = 5.9

beta2 = 25

thk2 = .156

B2 = .47

lamba2 = 1 - .1


#====================================================
#======== STEADY OPERATING STATE PARAMETERS =========
#====================================================

#Generated parameters entrance
#==============================

#metal bloackage
Bk1 = (thk*B1*z)/(math.sin(beta1*np.pi/180))

#
A1 = (np.pi*B1*(Rhub1 + Rtip1) - Bk1)*lamba1


#Velocity triangles
#=================================

#meridional velocity
Cm1 = (144*m)/(rho1*A1)

#blade tangnetial velocity
U = (2*np.pi*Rtip1*N)/(720)

Cu1 = Cm1/math.tan(alpha*np.pi/180)

C1 = (Cm1**2 + Cu1**2)**.5

Beta_f1 = math.atan(Cm1/(Cu1 - C1))*180/np.pi

#incidence angle
i1 = beta1 - Beta_f1

#Realtive Velocity triangles
#=================================

#realtive tangential velocity
Wu1 = Cm1*math.tan((90 - Beta_f1)*np.pi/180)

#relative velocity
W1 = (Wu1**2 + Cm1**2)**.5




delta = Rhub1/Rtip1 + Rhub2/Rtip2

X = S + 2*Rhub2/Rtip2


#slip factor
sigma = 1/(1 + ((1 + 0.6*math.sin(beta2*np.pi/180)) / (z*(1 + delta)*X**2 + 0.25*(1 - delta)**2)**.5))

#Generated parameters exit
#==============================

Bk2 = (thk2*B2*z)/(math.sin(beta2*np.pi/180))

A2 = (np.pi*B2*(Rtip2+Rhub2) - Bk2)*lamba2

Cm2 = (144*m)/(rho1*A2)

#Head rise (loop)
#============================
H2 = 900
H2i = 850
yd = H2/H2i
U_int = 5

U2 = U_int
x = 1
for x in range(800):
    U2new = U2 - ((((U2**(2) + U2*Cm2*math.tan(beta2*np.pi/180)+ U2**(2)*(1 - sigma) - U*Cu1)*yd/g)-H2)/((2*U2 + Cm2*math.tan(beta2*np.pi/180)+ 2*U2*(1 - sigma))*yd/g))
    U2 = U2new
    print(f"U2: {U2:.4f}")
 
#Generated parameters exit (P2)
#==============================
Wu2 = Cm2*math.tan(beta2*np.pi/180)+ U2*(1 - sigma)
    
    
Cu2 = U2 + Wu2
   
Beta_f2 = math.atan(Cm2/Wu2)


#rotor effciency calculations
#=============================================

#Specfic speed at BEP
Ns = (np.pi*N*Q**.5)/(30*g**.75 * H2**.75) 

nhyd_d = 0.41989 + 2.1524*Ns - 3.1434*Ns**2  + 1.5673*Ns**3


#pressure at exit point 2
Pt2 = (H2 * rho1)/144 + Pt1

NPSH = (144*(Pt1 - Pv))/(rho1)

Nss = (N*Q**.5)/(NPSH)**.75

#Diffusion System
#=================================

#loading parameter Cavitation checks
BB = 1.2
Ps_throat = Pt1 - ((C1*BB)**2 * rho1)/(2*144*g)

#total pressure loss coeffcient diffuser 0.2-0.9
C2 = (Cm2**2 + Cu2**2)**.5

Ps_2 = Pt2 - ((C2*BB)**2 *rho1)/(2*144*g)

W_loss = .5

Pt4 = -W_loss*(Pt2 - Ps_2) + Pt2

A_throat = A2

A4 = .1

Cthroat = (144*m)/(rho1*A_throat)

C4 = (144*m)/(rho1*A4)
#L = Cthroat/()

Ps4 = Pt4 - (C4**2 * rho1)/(2*g)

H4 = ((Pt4 - Pt1)*144)/rho1


#disk pumpage windage
Hpd = 32* K* N**(3) * (Rhub2/12)**(5)


#volume flow effciency
n_vol = m/(m+mL)

n_mech = .98

HP = (m*H2)/(550*(n_vol*yd*n_mech)) + Hpd

#Overall Pump effciency
n_4 = (m*H4)/(550*HP)




#Off Performance Design Parameters
#=======================================


#Off design Suction perfroamnce
#=====================================
f = 0
nss_l = []
f_l = []
for i in range (15):
    f = 0 + .1*i
    nssr = (-0.28607 + 4.14245*f - 12.0967*f**2 + 20.708*f**3 - 15.42122*f**4 + 3.9366*f**5)/Nss
    nss_l.insert(i, nssr)
    f_l.insert(i, f)
    

fig = plt.figure()
x = f_l
y = nss_l

plt.plot(x, y)
plt.show()


#Off design effieciency 
#=====================================
#rotor effciency

f = 0
nhyd = []
f_l = []
for i in range (15):
    f = 0 + .1*i
    nyhde =  (0.86387 + 0.3096*f - 0.14086*f**2 - 0.029265*f**3)
    nhyd.insert(i, nyhde)
    f_l.insert(i, f)



fig = plt.figure()
x = f_l
y = nhyd

plt.plot(x, y)
plt.show()



#Off design rotor head slip factor
#=====================================

f = 0
sigma_off = []
f_l = []
for i in range (15):
    f = 0 + .1*i
    sigma_r =   1.534988 - 0.6681668*f + 0.077472*f**2 + 0.0571508*f**3
    sigma_off.insert(i, sigma_r)
    f_l.insert(i, f)



fig = plt.figure()
x = f_l
y = sigma_off

plt.plot(x, y)
plt.show()



#Off design pressure loss coeffcient
#=====================================

f = 0
W_off = []
L = []
for i in range (15):
    l = 0 + .1*i
    W =   1.8151 - 1.83527*l + 0.8798*l**2 + 0.18765*l**3
    W_off.insert(i, W)
    L.insert(i, l)



fig = plt.figure()
x = f_l
y = W_off

plt.plot(x, y)
plt.show()
