import numpy as np
import constants as const
import equation_of_state
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from derivatives import derivs

def shootf(trial_params, M_star, X, Y):
    """
    trial_params: [Pc, Tc, R_star, L_star] [cite: 24]
    This function calculates the mismatch at the fitting point. 
    """
    Pc, Tc, R_star, L_star = trial_params
    m_fit = 0.5 * M_star # Common choice for fitting point
    
    # 1. Outward Integration (load1) [cite: 21, 22]
    # Start at m_epsilon (close to 0) to avoid singularities 
    y_center = load1(Pc, Tc, m_epsilon) 
    outward_sol = solve_ivp(derivs, [m_epsilon, m_fit], y_center, args=(X, Y,))
    
    # 2. Inward Integration (load2) [cite: 21, 23]
    y_surface = load2(R_star, L_star, M_star)
    inward_sol = solve_ivp(derivs, [M_star, m_fit], y_surface, args=(X, Y,))
    
    # 3. Mismatch 
    # Difference between outward_sol.y[:,-1] and inward_sol.y[:,-1]
    return outward_sol.y[:, -1] - inward_sol.y[:, -1]