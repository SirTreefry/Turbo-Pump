# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 20:53:49 2026

@author: Andrew Trefry
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

m = 5

mL = .1

rho1 = 62.4

N = 36000

Q = m/rho1

Pt1 = 45

Pv = .3

K = 3*10**(-12)

#Entrance or stage 1 values
#=======================================

#blade span hub to tip


#blade angle
beta1 = 25

alpha = 22

Rhub1 = .45

Rtip1 = .75

#blade span hub to tip (only calcuklated this way at entrance of impeller)
B1 = Rtip1 - Rhub1

lamba1 = 1 - .02

S = .1

#Entrance or stage 2 values
#=======================================
Rhub2 = .9

Rtip2 = .9

beta2 = -25

thk2 = .156

B2 = .35

lamba2 =  1 - .1

#Entrance or stage 3 values
#=======================================
A3 = 1 #*A3

A4 = 2.5

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

Beta_f1 = math.atan2(Cm1, U - Cu1) * 180 / np.pi

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

Bk2 = (thk2*B2*z)/(math.sin(abs(beta2)*np.pi/180))

A2 = (np.pi*B2*(Rtip2+Rhub2) - Bk2)*lamba2

Cm2 = (144*m)/(rho1*A2)


yd = .8
U_int = 5

U2 = (2*np.pi*Rtip2*N)/720

#Generated parameters exit (P2)
#==============================
Wu2 = Cm2*math.tan(beta2*np.pi/180)+ U2*(1 - sigma)
    
    
Cu2 = U2 + Wu2
   
Beta_f2 = math.atan2(Cm2,Wu2)


H2 = (U2*Cu2 - U*Cu1)*yd/g

#rotor effciency calculations
#=============================================

#Specfic speed at BEP
Ns = 21.2*N*(Q)**.5 /(H2)**.75




#pressure at exit point 2
Pt2 = (H2 * rho1)/144 + Pt1

NPSH = (144*(Pt1 - Pv))/(rho1)

Nss = (N*Q**.5)/(NPSH)**.75

#Diffusion System
#=================================
#NEEDS VERIFICATION IDENTRIFICATION CHECK REAL DESIGN
R3 = 1.25
A3 = A2

Cm3 = Cm2

Cu3 = Cu2*Rtip2/R3

#throat systems
#==========================================
A_throat = A2*.07
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


#disk pumpage windage
Hpd = 32* K* N**(3) * (Rhub2/12)**(5)


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

N_list = [10000, 15000, 20000, 25000, 30000, 36000, 40000]

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

        delta = Rhub1/Rtip1 + Rhub2/Rtip2


        X = S + 2*Rhub2/Rtip2


        #slip factor
        sigma = (1/(1 + ((1 + 0.6*math.sin(abs(beta2)*np.pi/180)) / (z*(1 + delta)*X**2 + 0.25*(1 - delta)**2)**.5)))

        #Generated parameters exit
        #==============================

        Bk2 = (thk2*B2*z)/(math.sin(abs(beta2)*np.pi/180))

        A2 = (np.pi*B2*(Rtip2+Rhub2) - Bk2)*lamba2

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
        

        L_design = Cthroat/(Cu3**2 + Cm3**2)**.5 
        print(L_design)
        w_off = 1.8151 - 1.83527*Ldesign + 0.8798*L_design**2 + 0.18765*L_design**3
        
        Pt4 = -W_loss*w_off*(Pt2 - Ps_2)+Pt2
        H4 = (Pt4 - Pt1)*144/g
        H_off.insert(j,H4)
        
        C4 = (144 * m) / (rho1 * A4) 

        # Convert Total Pressure (Pt4) to Static Pressure (Ps4)
        Ps4 = Pt4 - (C4**2 * rho1) / (2 * 144 * g)
        
        P_off.insert(j,Pt4)
        if Ps_2 < 0:
            print("FAILED at")
            print("Q =",Q)
            print("RPM:",N_list[k])
            print("Pt2 =", Pt2)
            print("C2  =", C2)
            print("Ps_2 =", Ps_2)
            
            
        if q >= Pt2:
            Ps_2 = 0.0
            # mark as non-physical / stalled
            break
        
        print("C1", C1)
        print("Pt1 =",Pt1)
        
        print("C2", C2)
        print("Pt2 =",Pt2)
        
        print("C4", C4)
        print("Pt4 =",Pt4)
        
        
        print("RPM:",N_list[k])
        print("Q", Q)
        

        
       
            


        
    
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

rpm_labels = ['RPM 10,000', 'RPM 15,000', 'RPM 20,000', 'RPM 30,000', 'RPM 25,000', 'RPM 36,000', 'RPM 40,000']

plt.figure(figsize=(8,6))

for Q, H, label in zip(Q_list, H_list, rpm_labels):
    Q = np.array(Q)
    H = np.array(H)
    # List to hold stall points
    stall_Q = []
    stall_H = []
    # Detach at j=8
    Q_detach = Q[10:]
    H_detach = H[10:]
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

rpm_labels = ['RPM 10,000', 'RPM 15,000', 'RPM 20,000', 'RPM 20,000', 'RPM 25,000', 'RPM 36,000', 'RPM 40,000']

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

rpm_labels = ['RPM 10,000', 'RPM 15,000', 'RPM 20,000', 'RPM 20,000', 'RPM 25,000', 'RPM 36,000', 'RPM 40,000']

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





