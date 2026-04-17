import numpy as np

def pp_chain(psi, f11, g11, rho, X, T9):
    # Calculate the energy generation rate for the pp chain
    epsilon_pp = 2.57e4*psi*f11*g11*rho*(X**2)*T9**(-2/3)*np.exp(-3.381/(T9**(1/3)))
    return epsilon_pp

def cno_cycle(g14_1, X_CNO, X, rho, T9):
    # Calculate the energy generation rate for the CNO cycle
    epsilon_CNO = 8.24*1e25*g14_1*X_CNO*X*rho*T9**(-2/3)*np.exp(-15.231*T9**(-1/3) - (T9/0.8)**2)
    return epsilon_CNO

def total_energy_generation(epsilon_pp, epsilon_CNO):
    # Calculate the total energy generation rate
    epsilon_total = epsilon_pp + epsilon_CNO
    return epsilon_total