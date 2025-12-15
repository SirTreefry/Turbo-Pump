# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 12:37:27 2025

@author: trefr
"""

#Inducer design for fuel pump
#===================================
import numpy as np

# -------- GIVEN DESIGN INPUTS --------
Q = 0.05            # volumetric flow rate [m^3/s]
N = 30000           # rotational speed [RPM]
rho = 1140          # fluid density [kg/m^3] (LOX example)
g = 9.81            # gravity [m/s^2]



#Blumfield criterion
#====================================
hub_tip_ratio = 0.5   # r_hub / r_tip

phi = 0.07          # flow coefficient (typical inducer range 0.05–0.10)(need to add eq later)
#for smaller than 0.1 look at - eq 63 for more accurate number can rearrnage for ss

K = (2*phi)/(1 - 2*phi) #blade tip cavation number

tau = K + K*phi**2 + phi**2  # conservative cavitation parameter (from SP-8052 trends)

Ss_correct = 8147*phi**.5 * tau**(-.75)
#suction specific speed
Ss = Ss_correct*(1 - hub_tip_ratio**2)**.5

q_correct = Q/(1 - hub_tip_ratio**2)

#Basic calculations
NPSH_a = 6.0          # required NPSH [m] 6-12 percent total head need to reclaulte based on head
NPSH_r = (N*q_correct/Ss_correct)**.75

D = 0.37843*(Q/((1 - hub_tip_ratio**2)*N*phi))**(1/3)

#Cavitation Properties
alpha = 

Z = 3                 # blade count (2–4 typical for inducers)


# Angular velocity
omega = 2 * np.pi * N / 60  # rad/s



U_tip = np.sqrt((2 * g * NPSH_a) / tau)
r_tip = U_tip / omega

#tau = 0.15   conservative cavitation parameter (from SP-8052 trends)

U_tip = np.sqrt((2 * g * NPSH_a) / tau)
r_tip = U_tip / omega


r_hub = hub_tip_ratio * r_tip

A = np.pi * (r_tip**2 - r_hub**2)


V_ax = Q / A

# Check flow coefficient consistency
phi_calc = V_ax / U_tip
print("Flow coefficient check:", phi_calc)


# Blade inlet angles (relative flow)
beta_tip = np.arctan(V_ax / U_tip)
beta_hub = np.arctan(V_ax / (omega * r_hub))

# Convert to degrees
beta_tip_deg = np.degrees(beta_tip)
beta_hub_deg = np.degrees(beta_hub)

eta_ind = 0.85   # optimistic early-stage inducer efficiency

delta_h = eta_ind * U_tip * V_ax / g
delta_p = rho * g * delta_h


circumference = 2 * np.pi * r_tip
pitch = circumference / Z

# Recommended chord ~ 0.7–1.2 pitch
chord = 0.9 * pitch


print("\n--- INDUCER GEOMETRY OUTPUT ---")
print(f"Tip radius [m]: {r_tip:.4f}")
print(f"Hub radius [m]: {r_hub:.4f}")
print(f"Axial velocity [m/s]: {V_ax:.2f}")
print(f"Tip blade angle [deg]: {beta_tip_deg:.2f}")
print(f"Hub blade angle [deg]: {beta_hub_deg:.2f}")
print(f"Blade count: {Z}")
print(f"Blade chord [m]: {chord:.4f}")
print(f"Estimated inducer pressure rise [Pa]: {delta_p:.0f}")


#Andrews Anylasis

#q_correct = 

#corrected_suction = 
