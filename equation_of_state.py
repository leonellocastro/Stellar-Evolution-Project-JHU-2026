import numpy as np
import constants
# Constants
a = constants.a # radiation constant in CGS
N_A = constants.N_A # Avogadro's number
k_B = constants.k_B
mu = constants.mu

# Equation of state in terms of P, T, X, Y

def equation_of_state(P,T):
    # P = P_gas + P_rad
    # P = (rho * k * T) / (mu(X,Y) * m_H) + (1/3) * a * T**4
    # Approximate mean molecular weight for a fully ionized gas
    rho = (P - (1/3) * a * T**4) * mu / (N_A * k_B * T)
    return rho