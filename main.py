# main.py

import numpy as np
from scipy.optimize import root
from scipy.optimize import least_squares
import pickle
import constants
from load_table import load_table
from find_exact_table_id import find_exact_table_id
from shootf import shootf, load1, load2
from derivs import derivs

delta_m = constants.delta_m
filename = constants.filename
M_target = constants.M_target

LSUN = constants.LSUN
RSUN = constants.RSUN

X = constants.X
Y = constants.Y
Z = constants.Z

# Initial Guesses

Pc_guess = constants.Pc_guess
Tc_guess = constants.Tc_guess
L_guess = constants.L_guess
R_guess = constants.R_guess

# Load opacity table and calculate average kappa
table_id = find_exact_table_id(filename, X, Y, Z)

table_data = load_table(filename, table_id)
kappa_matrix = table_data[2]

log_kappa_avg = np.nanmean(kappa_matrix) 
kappa_avg = 10**log_kappa_avg

start_out = load1(delta_m, Pc_guess, Tc_guess)
start_in = load2(M_target, L_guess, R_guess, kappa_avg)

# Example
M_target = constants.M_target
delta_m = constants.delta_m

# Opacity Matrix

# 1. Fetch the matrix for the specific composition
# load_table returns (lT, lR, kappa_mat, X, Y, Z)
table_data = load_table(filename, table_id)
kappa_matrix = table_data[2]

# 2. Calculate the average kappa
log_kappa_avg = np.nanmean(kappa_matrix) 
kappa_avg = 10**log_kappa_avg

start_out = load1(delta_m, Pc_guess, Tc_guess)
start_in = load2(M_target, L_guess, R_guess, kappa_avg)

# Calculate Initial Conditions at Both Boundaries
print("Outward Start (L, P, r, T):", start_out)
print("Inward Start  (L, P, r, T):", start_in)

# Inner Boundary (Center)

derivs_inner = derivs(delta_m, start_out, X, Y, Z, filename)

# Outer Boundary (Surface)

derivs_outer = derivs(M_target, start_in, X, Y, Z, filename)

print(f"{'Variable':<12} | {'Inner Deriv':<20} | {'Outer Deriv':<20}")
print("-" * 60)

vars = ['dl/dm (erg/s/g)', 'dp/dm (dyn/g)  ', 'dr/dm (cm/g)   ', 'dt/dm (K/g)    ']
for i in range(4):
    print(f"{vars[i]:<12} | {derivs_inner[i]:<20.4e} | {derivs_outer[i]:<20.4e}")

# Initial guesses for [Pc, Tc, L_star, R_star]
initial_guesses = [Pc_guess, Tc_guess, L_guess, R_guess]

# Use 'root' to find where shootf returns zero
sol = root(shootf, initial_guesses, args=(M_target, delta_m, X, Y, Z, kappa_avg, filename), method='hybr')

"""
# Define bounds: (min_Pc, min_Tc, min_L, min_R), (max_Pc, max_Tc, max_L, max_R)
# Adjust these based on 5 Msun expectations
lower_bounds = [1e15, 1e6, 1e32, 5e10]
upper_bounds = [1e18, 5e7, 1e36, 5e11]

res = least_squares(shootf, initial_guesses, args=(M_target, delta_m, X, Y, Z, kappa_avg, filename), bounds=(lower_bounds, upper_bounds), ftol=1e-7)
sol = res.x
"""

if sol.success:
    final_Pc, final_Tc, final_L, final_R = sol.x
    print(f"Converged! L={final_L/LSUN:.2f} L_sun, R={final_R/RSUN:.2f} R_sun, Pc={final_Pc:.2e} dyn/cm^2, Tc={final_Tc:.2e} K")
else:
    print("Convergence failed. Check initial guesses.")

if sol.success:
    save_data = {'sol': sol, 'kappa_avg': kappa_avg}
    with open('model_result.pkl', 'wb') as f:
        pickle.dump(save_data, f)
    print("Solution and kappa_avg saved to model_result.pkl")