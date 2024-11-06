import numpy as np
import matplotlib.pyplot as plt

g0_f_H2O = -237180  # J/mol for H2O at 298K and 1 atm (given at the end)
g0_f_H2 = -100000 # J/mol for H2 (recheck)
g0_f_O2 = 0  # J/mol for O2 (recheck)

n_electrons = 2 # For hydrogen oxidation
F = 96485 #given data
R = 8.314  # Universal gas constant (standard value)
T_stack = 353 # Stack temperature (on matlab parameters) but given 80C(353K) on the intro but you get negative curve for that temperature

ne_a = 4 #number of elec in anode
ne_c =1
beta = 0.5 #symmetry factor(given)

k_a = 1 
k_c = 1
i_loss= 0.002 


C1, C2 = 180, 16.4 #given in chapter 3

active_area = 232 

def gibbs_free_energy(g0_f_H2O, g0_f_H2, g0_f_O2):
    delta_g_f = g0_f_H2O - g0_f_H2 - 0.5 * g0_f_O2
    return delta_g_f

def theoretical_cell_potential(delta_g_f):
    return -delta_g_f / (n_electrons * F)

def nernst_potential(E0, T=T_stack, a_H2=1, a_O2=0.21, a_H2O=1): #so the value of H20 is given but H2 and O2 is take reference from dipesh's code , need to verify
    return E0 + (R * T / (2 * F)) * np.log(a_H2 * np.sqrt(a_O2) / a_H2O)

#Gibbs free energy and theoretical cell potential first page equation and formulas
delta_g_f = gibbs_free_energy(g0_f_H2O, g0_f_H2, g0_f_O2)
E0 = theoretical_cell_potential(delta_g_f)

# Activation losses (Tafel equation), formula used from intro part 
def activation_losses(current_density, i0_a=0.1, i0_c=0.1, alpha_a=2, alpha_c=0.5):
    act_loss_anode = (R * T_stack / (alpha_a * F)) * np.log((current_density + i_loss) / i0_a )
    act_loss_cathode = (R * T_stack / (alpha_c * F)) * np.log((current_density + i_loss) / i0_c)
    return act_loss_anode + act_loss_cathode

#Ohmic loss  (formula is in page 14 chapter 3)
def ohmic_losses(current_density, lambda_mem, t_mem=0.1):
    resistance = C1 * (1 + 0.03 * current_density + 0.062 * T_stack / 303 * current_density ** 2.5) / (
        lambda_mem - 0.634 - 3 * current_density) * np.exp(C2 * (T_stack - 303) / T_stack) * t_mem
    return resistance * current_density

def cell_voltage(E0, current_density, lambda_mem):
    nernst = nernst_potential(E0)
    act_losses = activation_losses(current_density)
    ohm_losses = ohmic_losses(current_density, lambda_mem)
    return nernst - act_losses - ohm_losses


alpha_a = beta * ne_a
alpha_c = (1-beta) * ne_c

i0_a = n_electrons * F * k_a * np.exp((1 - beta) * n_electrons * F * E0 / (R * T_stack))
i0_c = n_electrons * F * k_c * np.exp(-beta * n_electrons * F * E0 / (R * T_stack))

current_densities = np.linspace(0.01, 1.0, 100)  # A/cm^2
lambda_mem = 14  #since the activity of water is given 1 

# Calculate voltages and corresponding currents
voltages = [cell_voltage(E0, i, lambda_mem) for i in current_densities]
currents = [i * active_area for i in current_densities]  # current density times the area (need to confirm the values but not needed if plotting V vs c_density)

# Plot the VI curve
plt.figure(figsize=(8, 6))
plt.plot(current_densities, voltages, label='PEM Fuel Cell VI Curve')
plt.xlabel("Current Desnsity (A/cm^2)")
plt.ylabel("Cell Voltage (V)")
plt.title("PEM Fuel Cell VI Curve")
plt.legend()
plt.grid(True)
plt.show()

print(E0)







