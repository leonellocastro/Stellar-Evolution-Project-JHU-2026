import numpy as np
import constants as const
import equation_of_state
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from derivatives import derivs

# STEP 1 (Mass and Composition)

# Define Mass of Star

M_star = 5*const.M_sun

# Define Composition of Star

X = 0.70
Y = 0.28
Z = 1 - X - Y

# 4. Global Driver 
# Use fsolve to find Pc, Tc, R_star, L_star that make shootf return zeros
initial_guesses = [log_Pc_guess, log_Tc_guess, R_guess, L_guess]
solution = fsolve(shootf, initial_guesses, args=(M_star, X, Y))