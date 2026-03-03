# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 00:17:21 2026

@author: trefr
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import sympy as sp
from sympy import cos, nsolve, Symbol
#1D Pump Sizing for impellers
#Based on Pumpas method
g = 32.2

#Inputs dimensions
#===========================================

#number of blades
z = 5


m = 15.309860702

mL = .1*m

rho1 = 64.2

N = 20000

Q = m/rho1 


Pt1 = 45

Pv = .3

K = 3*10**(-12)

#Entrance or stage 1 values
#=======================================

#blade span hub to tip


#blade angle
beta1 = 21.9

alpha = 90 #should be zero

Rhub1 = .2

Rtip1 = .6

#blade span hub to tip (only calcuklated this way at entrance of impeller)
B1 = Rtip1 - Rhub1

lamba1 = 1 - .2

#normal blade thickeness
thk = .008

S = .5

#Entrance or stage 2 values
#=======================================


beta2 = 20

thk2 = .008

B2 = .15

lamba2 =  1 - .25

#Entrance or stage 3 values
#=======================================
A3 = 1 #*A3

A4 = 2

#====================================================
#======== STEADY OPERATING STATE PARAMETERS =========
#====================================================

#Generated parameters entrance
#==============================

#metal bloackage
Bk1 = (thk*B1*z)/(math.sin(abs(beta1)*np.pi/180))

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

Beta_f1 = math.atan2(Cm1, U - Cu1) * 180 / np.pi

#incidence angle
i1 = beta1 - Beta_f1

#Realtive Velocity triangles
#=================================

#realtive tangential velocity
Wu1 = Cm1*math.tan((90 - Beta_f1)*np.pi/180)

#relative velocity
W1 = (Wu1**2 + Cm1**2)**.5




delta = Rhub1/Rtip1 + 1

X = S + 2


#slip factor
sigma = 1/(1 + ((1 + 0.6*math.sin(beta2*np.pi/180)) / (z*(1 + delta)*X**2 + 0.25*(1 - delta)**2)**.5))

#Generated parameters exit
#==============================

Bk2 = (thk2*B2*z)/(math.sin(abs(beta2)*np.pi/180))


import scipy.optimize as opt

def W1_sq(Rtip1):
    A1 = np.pi*(Rtip1**2 - Rhub1**2)  # simplified, no blockage
    Cm1 = (144*m)/(rho1*A1)
    U1 = (2*np.pi*Rtip1*N)/720
    return Cm1**2 + U1**2

result = opt.minimize_scalar(W1_sq, bounds=(0.3, 2.0), method='bounded')
Rtip1_opt = result.x


#rotor effciency calculations
#=============================================
H2 = 1653
#Specfic speed at BEP
Ns = (np.pi*(N)*(Q)**.5) /(30* g**.75  *H2**.75)

NPSH = (144*(Pt1 - Pv))/(rho1)

Nss = (N*Q**.5)/(NPSH)**.75


#Velocity U2 (loop)
#============================
yd =  0.41989 + 2.1524*Ns - 3.1434*Ns**2 + 1.5673 *Ns**3

#-(144*m*rho1*np.pi*B2*2*lamba2)/(rho1*(np.pi*B2*(2*Rtip2) - Bk2)*lamba2)**2
#((2*np.pi*Rtip2*N)/720)*((144*m)/(rho1*(np.pi*B2*(2*Rtip2) - Bk2)*lamba2))


# ((4*np.pi*N)/720)*((2*np.pi*Rtip2*N)/720)*((144*m)/(rho1*(np.pi*B2*(2*Rtip2) - Bk2)*lamba2))    + ((2*np.pi*Rtip2*N)/720)*-(144*m*rho1*np.pi*B2*2*lamba2)/(rho1*(np.pi*B2*(2*Rtip2) - Bk2)*lamba2)**2
Rtip2 = 0


def tip_radius():
    
    Rtip2 = Symbol('Rtip2')
    U2 = (2 * np.pi * Rtip2 * N) / 720
     

    #=============== euler values =============================

    Cm2 = (144*m)/(rho1*(np.pi*B2*(Rtip2*2) - Bk2)*lamba2)
    Cm2d = -(144*m*rho1*np.pi*B2*2*lamba2)/(rho1*(np.pi*B2*(2*Rtip2) - Bk2)*lamba2)**2
    Rtip = 8
    
    #================= slip factor ===========================
    Rrms2 = ((Rtip2**2 + Rtip2**2)/2)**0.5
    Rrms1 = ((Rtip1**2 + Rhub1**2)/2)**0.5
    delta = Rrms1 + Rrms2

    X = S + 2*Rrms2


    sigma = 1/(1 + ((1 + 0.6*sp.sin(beta2*np.pi/180)) / (z*(1 + delta)*X**2 + 0.25*(1 - delta)**2)**.5))
    
    #============== differntation bottom euler head ===========
    
    
    #f = ((U2**(2) + U2*Cm2*sp.tan(beta2*np.pi/180)+ U2**(2)*(1 - sigma) - U*Cu1)*yd/g)-H2
    f = ((U2**(2) - U2*Cm2*sp.tan(beta2*np.pi/180)- U2**(2)*(1 - sigma) - U*Cu1)*yd/g)-H2
    df = sp.diff(f,Rtip2)
    
    for x in range(100):
        
        
        
        
        
        
       
        fval = f.subs(Rtip2,Rtip)
        df_sub = df.subs(Rtip2,Rtip)
        
        #Rnew = Rtip - (((U2**(2) + U2*Cm2*sp.tan(beta2*np.pi/180)+ U2**(2)*(1 - sigma) - U*Cu1)*yd/g)-H2)/(df_sub)
        Rnew = Rtip - (fval)/(df_sub)
        
        sigma_val = float(sigma.subs(Rtip2, Rtip))
        print(sigma_val)
        
        Rtip = Rnew.subs(Rtip2,Rtip)
        
        
    
    Rtip2 = Rtip
    print(Rtip2)
    
    #old func what de hellli
    #Rnew = ((2*np.pi*Rtip2*N)/720) - ((( ((2*np.pi*Rtip2*N)/720)**(2) + ((2*np.pi*Rtip2*N)/720)*((144*m)/(rho1*(np.pi*B2*(2*Rtip2) - Bk2)*lamba2))*math.tan(beta2*np.pi/180)+ ((2*np.pi*Rtip2*N)/720)**(2)*(1 - sigma) - U*Cu1)*yd/g)-H2)/ ((2*((2*np.pi*N)/720) + ((4*np.pi*N)/720)*((2*np.pi*Rtip2*N)/720)*((144*m)/(rho1*(np.pi*B2*(2*Rtip2) - Bk2)*lamba2))*math.tan(beta2*np.pi/180)    + ((2*np.pi*Rtip2*N)/720)*(-144*m*rho1*np.pi*B2*2*lamba2)/(rho1*(np.pi*B2*(2*Rtip2) - Bk2)*lamba2)**2 *math.tan(beta2*np.pi/180)+ (2*((2*np.pi*N)/720)*(1 - sigma)))*yd/g)
    return Rtip2,sigma_val

Rtip2, sigma = tip_radius()

    
A2 = (np.pi*B2*(Rtip2*2) - Bk2)*lamba2
U2 = (2*np.pi*Rtip2*N)/720
print("u1",U)
print("cu1",Cu1)


#pressure at exit point 2
Pt2 = (H2 * rho1)/144 + Pt1

#Generated parameters exit (P2)
#==============================
Cm2 = (144*m)/(rho1*A2)
Wu2 = Cm2*math.tan(beta2*np.pi/180)+ U2*(1 - sigma)
 
    
Cu2 = U2 + Wu2
   
Beta_f2 = math.atan2(Cm2,Wu2)


Head_coeff = (g*H2)/(U2**2)

Flow_coeff = (Q)/(U2*(A2/144))

print("U2",U2)

print("sigma",sigma)
print("CU2", Cu2)
print("Beta2",beta2)
print("wu2 is",Wu2)
print("Wu1", Wu1)

print("sigma:",sigma)
print(f"A1: {A1:.4f}")
print(f"A2: {A2:.4f}")
print(f"Bk1: {Bk1:.4f}")
print(f"Bk2: {Bk2:.4f}")
print(f"Cm1: {Cm1:.4f}")
print(f"Cm2: {Cm2:.4f}")
print(f"Q: {Q:.4f}")
print("Effciency is :",yd)
print("Ns is :",Ns)
print("Nss is :",Nss)
print("euler Head Coeff is:",Head_coeff)
print("Flow coeff is:",Flow_coeff)
#Diffusion System
#=================================
#NEEDS VERIFICATION IDENTRIFICATION CHECK REAL DESIGN
R3 = 2


Cm3 = Cm2

Cu3 = Cu2*Rtip2/R3

#throat systems
#==========================================
A_throat = A2*.1
phi = Cm1/Q
Cthroat = (144*m)/(rho1*A_throat)

Ldesign = Cthroat/(Cu3**2 + Cm3**2)**.5 

#loading parameter Cavitation checks
BB = 1.2
Ps_throat = Pt1 - ((C1*BB)**2 * rho1)/(2*144*g)

#total pressure loss coeffcient diffuser 0.2-0.9
C2 = (Cm2**2 + Cu2**2)**.5

Ps_2 = Pt2 - ((C2)**2 *rho1)/(2*144*g)

W_loss = .15

Pt4 = -W_loss*(Pt2 - Ps_2) + Pt2




C4 = (144*m)/(rho1*A4)
#L = Cthroat/()

Ps4 = Pt4 - (C4**2 * rho1)/(2*g)

H4 = ((Pt4 - Pt1)*144)/rho1


Head_coeff = (g*H4)/(U2**2)
print("Head Coeff is:",Head_coeff)
print("total outlet pressure",Pt4)
#disk pumpage windage
Hpd = 32* K* N**(3) * (Rtip2/12)**(5)


#volume flow effciency
n_vol = m/(m+mL)

n_mech = .98

HP = (m*H2)/(550*(n_vol*yd*n_mech)) + Hpd

#Overall Pump effciency
n_4 = (m*H4)/(550*HP)





#====================================================
#======== Off Performance Design Parameters =========
#====================================================


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


#Pump Curve Headrise verus flow rate
#===============================================

N_list = [10000, 15000,20000, 22000, 25000,27000,30000]

#Nested loops for total graoh changed every RPM loop
H_tot = []
Q_tot = []
P_tot = []


#loop values re-placed every perfrmmance iteration
H_off = []
f = 0
Q_off = []
m_off = []
P_off = []
L_off = []
Q_base = Q

for k in range (7):
    print(N_list[k]/N)
    Q_it = Q_base*(N_list[k]/N)
   
    Q_off = [] 
    H_off = []
    m_off = []
    P_off = []
    L_off = []
    j = 0
    for j in range (15):
        
        i = 1 + j
        
        f = 0 + .1*i
        
        Q = Q_it*f
        
        
    
        Q_off.insert(j, Q)
    
        m = Q*rho1
    
        m_off.insert(j,m)
    
    
    
        #meridional velocity
        Cm1 = (144*m)/(rho1*A1)

        #blade tangnetial velocity
        U = (2*np.pi*Rtip1*N_list[k])/(720)

        Cu1 = Cm1/math.tan(alpha*np.pi/180)

        C1 = (Cm1**2 + Cu1**2)**.5

        Beta_f1 =  math.atan2(Cm1, U - Cu1) * 180 / np.pi

        #incidence angle
        i1 = beta1 - Beta_f1

        #Realtive Velocity triangles
        #=================================

        #realtive tangential velocity
        Wu1 = Cm1*math.tan((90 - Beta_f1)*np.pi/180)

        #relative velocity
        W1 = (Wu1**2 + Cm1**2)**.5

        delta = Rhub1/Rtip1 + Rtip1/Rtip2


        X = S + 2*Rtip2/Rtip2


        #slip factor
        sigma = (1/(1 + ((1 + 0.6*math.sin(abs(beta2)*np.pi/180)) / (z*(1 + delta)*X**2 + 0.25*(1 - delta)**2)**.5)))

        #Generated parameters exit
        #==============================

        Bk2 = (thk2*B2*z)/(math.sin(abs(beta2)*np.pi/180))

        A2 = (np.pi*B2*(Rtip2+Rtip2) - Bk2)*lamba2

        Cm2 = (144*m)/(rho1*A2)
        
        U2 = (2*np.pi*Rtip2*N_list[k])/720
        
        Wu2 = Cm2*math.tan(beta2*np.pi/180) + U2*(1 - sigma*sigma_off[j]) #
       
        Cu2 = U2 + Wu2
        H = (U2*Cu2 - U*Cu1)*(yd*nhyd[j])/g
        
        Pt2 = (H * rho1)/144 + Pt1
        
        
        
        
        
        
        
        
      

    
        
        C2 = (Cu2**2 + Cm2**2)**.5
        Ps_2 = Pt2 - (C2**2 *rho1)/(2*144*g)
       
     
        q = (rho1 * C2**2) / (2 * 144 * g)
        
        Cthroat = 144*m/(rho1*A_throat)
        
        Cm3 = Cm2

        Cu3 = Cu2*Rtip2/R3
        L_design = Cthroat/(Cu3**2 + Cm3**2)**.5 
        print(L_design)
        w_off = 1.8151 - 1.83527*L_design + 0.8798*L_design**2 + 0.18765*L_design**3
        
        Pt4 = -W_loss*w_off*(Pt2 - Ps_2)+Pt2
        H4 = (Pt4 - Pt1)*144/rho1
        H_off.insert(j,H4)
        
        C4 = (144 * m) / (rho1 * A4) 

        # Convert Total Pressure (Pt4) to Static Pressure (Ps4)
        Ps4 = Pt4 - (C4**2 * rho1) / (2 * 144 * g)
        
        P_off.insert(j,Pt4)
        if Ps_2 < 0:
            print()
            
            
            
        if q >= Pt2:
            Ps_2 = 0.0
            # mark as non-physical / stalled
            break
        
        #print("C1", C1)
        print(f"F={f:.2f}, Q={Q:.4f}, H_rotor={H:.1f}, Pt2={Pt2:.2f}, Ps2={Ps_2:.2f}, L={L_design:.3f}, Pt4={Pt4:.2f}, H4={H4:.1f}")
        

        
       
            


        
    
    Q_tot.insert(k,Q_off)
    H_tot.insert(k,H_off)
    P_tot.insert(k,P_off)
    
Q_0 = Q_tot[0]
Q_1 = Q_tot[1]
Q_2 = Q_tot[2]
Q_3 = Q_tot[3]
Q_4 = Q_tot[4]
Q_5 = Q_tot[5]
Q_6 = Q_tot[6]

H_0 = H_tot[0]
H_1 = H_tot[1]
H_2 = H_tot[2]
H_3 = H_tot[3]
H_4 = H_tot[4]
H_5 = H_tot[5]
H_6 = H_tot[6]

P_0 = P_tot[0]
P_1 = P_tot[1]
P_2 = P_tot[2]
P_3 = P_tot[3]
P_4 = P_tot[4]
P_5 = P_tot[5]
P_6 = P_tot[6]


plt.plot(Q_0, H_0, label='RPM 10,000')
plt.plot(Q_1, H_1, label='RPM 15,000')
plt.plot(Q_2, H_2, label='RPM 20,000')
plt.plot(Q_3, H_3, label='RPM 20,000')
plt.plot(Q_4, H_4, label='RPM 25,000')
plt.plot(Q_5, H_5, label='RPM 36,000')
plt.plot(Q_6, H_6, label='RPM 40,000')




plt.xlabel('Flow Q (Ib/s)')
plt.ylabel('Head (ft)')
plt.title('Mass Flow Evolution in Reactors')
plt.legend()
plt.grid(True)
plt.show()    




# Your existing Q and H lists
Q_list = [Q_0, Q_1, Q_2, Q_3, Q_4, Q_5, Q_6]
H_list = [H_0, H_1, H_2, H_3, H_4, H_5, H_6]

rpm_labels = ['RPM 5,000', 'RPM 10,000', 'RPM 15,000', 'RPM 20,000', 'RPM 24,500', 'RPM 27,000', 'RPM 30,000']

plt.figure(figsize=(8,6))

for Q, H, label in zip(Q_list, H_list, rpm_labels):
    Q = np.array(Q)
    H = np.array(H)
    # List to hold stall points
    stall_Q = []
    stall_H = []
    # Detach at j=8
    Q_detach = Q[1:]
    H_detach = H[1:]
    # Add the first point to the stall line
    stall_Q.append(Q_detach[0])
    stall_H.append(H_detach[0])
    
    # Plot the detached curve
    plt.plot(Q_detach, H_detach, label=label)

# Plot the stall line
plt.plot(stall_Q, stall_H, 'k--', linewidth=2, label='Stall Line')  # black dashed line

plt.xlabel('Flow Q (lb/s)')
plt.ylabel('Head (ft)')
plt.title('Pump Map with Stall Line')
plt.legend()
plt.grid(True)
plt.show()



# Your existing Q and H lists
Q_list = [Q_0, Q_1, Q_2, Q_3, Q_4, Q_5, Q_6]
H_list = [H_0, H_1, H_2, H_3, H_4, H_5, H_6]

rpm_labels = ['RPM 5,000', 'RPM 10,000', 'RPM 15,000', 'RPM 20,000', 'RPM 24,500', 'RPM 27,000', 'RPM 30,000']

plt.figure(figsize=(8,6))

stall_Q = []
stall_H = []

# loop over all curves
for Q, H, label in zip(Q_list, H_list, rpm_labels):
    Q = np.array(Q)
    H = np.array(H)
    
    # check if curve is long enough to detach at j=8
    if len(Q) > 8:
        Q_detach = Q[8:]
        H_detach = H[8:]
    else:
        Q_detach = Q
        H_detach = H
    
    # store the first point of the detached curve for the stall line
    stall_Q.append(Q_detach[0])
    stall_H.append(H_detach[0])
    
    # plot the detached curve
    plt.plot(Q_detach, H_detach, label=label)

# convert stall lists to arrays for sorting (optional, ensures line looks nice)
stall_Q = np.array(stall_Q)
stall_H = np.array(stall_H)

# sort by flow to make the stall line connect properly
sorted_indices = np.argsort(stall_Q)
stall_Q = stall_Q[sorted_indices]
stall_H = stall_H[sorted_indices]

# plot the stall line
plt.plot(stall_Q, stall_H, 'k--', linewidth=2, label='Stall Line')

plt.xlabel('Flow Q (lb/s)')
plt.ylabel('Head (ft)')
plt.title('Pump Map')
plt.legend()
plt.grid(True)
plt.show()





  


# Your existing Q and H lists
Q_list = [Q_0, Q_1, Q_2, Q_3, Q_4, Q_5, Q_6]
H_list = [P_0, P_1, P_2, P_3, P_4, P_5, P_6]

rpm_labels = ['RPM 5,000', 'RPM 10,000', 'RPM 15,000', 'RPM 20,000', 'RPM 24,500', 'RPM 27,000', 'RPM 30,000']

plt.figure(figsize=(8,6))

stall_Q = []
stall_H = []

# loop over all curves
for Q, H, label in zip(Q_list, H_list, rpm_labels):
    Q = np.array(Q)
    H = np.array(H)
    
    # check if curve is long enough to detach at j=8
    if len(Q) > 8:
        Q_detach = Q[8:]
        H_detach = H[8:]
    else:
        Q_detach = Q
        H_detach = H
    
    # store the first point of the detached curve for the stall line
    stall_Q.append(Q_detach[0])
    stall_H.append(H_detach[0])
    
    # plot the detached curve
    plt.plot(Q_detach, H_detach, label=label)

# convert stall lists to arrays for sorting (optional, ensures line looks nice)
stall_Q = np.array(stall_Q)
stall_H = np.array(stall_H)

# sort by flow to make the stall line connect properly
sorted_indices = np.argsort(stall_Q)
stall_Q = stall_Q[sorted_indices]
stall_H = stall_H[sorted_indices]

# plot the stall line
plt.plot(stall_Q, stall_H, 'k--', linewidth=2, label='Stall Line')

plt.xlabel('Flow Q (lb/s)')
plt.ylabel('Pressure (Psia)')
plt.title('Pump Map')
plt.legend()
plt.grid(True)
plt.show()







rpm_labels = ['RPM 5,000', 'RPM 10,000', 'RPM 15,000', 'RPM 20,000', 'RPM 24,500', 'RPM 27,000', 'RPM 30,000']

# ---- HEAD MAP ----
fig, ax = plt.subplots(figsize=(9, 6))
stall_Q, stall_H = [], []

for Q_c, H_c, lbl in zip(Q_list, H_list, rpm_labels):
    Q_c = np.array(Q_c)
    H_c = np.array(H_c)
    if len(H_c) == 0:
        continue
    
    # find BEP (max head) and plot from start to 20% past BEP
    bep_idx = int(np.argmax(H_c))
    end_idx = bep_idx
    for i in range(bep_idx, len(H_c)):
        if H_c[i] < 0.80 * H_c[bep_idx]:
            end_idx = i
            break
        end_idx = i
    
    ax.plot(Q_c[bep_idx:end_idx+1], H_c[bep_idx:end_idx+1], linewidth=2, label=lbl)
    ax.scatter(Q_c[bep_idx], H_c[bep_idx], marker='*', s=150, color='black', zorder=5)
    stall_Q.append(Q_c[1])
    stall_H.append(H_c[1])

# stall line through leftmost points
arr = sorted(zip(stall_Q, stall_H))
sq, sh = zip(*arr)
ax.plot(sq, sh, 'k--', linewidth=2, label='Stall Line')
ax.scatter([], [], marker='*', s=150, color='black', label='BEP')
ax.set_xlabel('Flow Q (lb/s)', fontsize=12)
ax.set_ylabel('Head (ft)', fontsize=12)
ax.set_title('Pump Map', fontsize=13)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()

# ---- PRESSURE MAP ----
Q_list_p = [Q_0, Q_1, Q_2, Q_3, Q_4, Q_5, Q_6]
P_list   = [P_0, P_1, P_2, P_3, P_4, P_5, P_6]

fig, ax = plt.subplots(figsize=(9, 6))
stall_Q, stall_P = [], []

for Q_c, P_c, lbl in zip(Q_list_p, P_list, rpm_labels):
    Q_c = np.array(Q_c)
    P_c = np.array(P_c)
    if len(P_c) == 0:
        continue

    bep_idx = int(np.argmax(P_c))
    end_idx = bep_idx
    for i in range(bep_idx, len(P_c)):
        if P_c[i] < 0.80 * P_c[bep_idx]:
            end_idx = i
            break
        end_idx = i

    ax.plot(Q_c[bep_idx:end_idx+1], H_c[bep_idx:end_idx+1], linewidth=2, label=lbl)
    ax.scatter(Q_c[bep_idx], P_c[bep_idx], marker='*', s=150, color='black', zorder=5)
    stall_Q.append(Q_c[1])
    stall_P.append(P_c[1])

arr = sorted(zip(stall_Q, stall_P))
sq, sp_ = zip(*arr)
ax.plot(sq, sp_, 'k--', linewidth=2, label='Stall Line')
ax.scatter([], [], marker='*', s=150, color='black', label='BEP')
ax.set_xlabel('Flow Q (lb/s)', fontsize=12)
ax.set_ylabel('Pressure (psia)', fontsize=12)
ax.set_title('Pump Map', fontsize=13)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()
