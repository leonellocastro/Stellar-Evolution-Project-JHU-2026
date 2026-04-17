# derivs.py (ODE derivatives function)

import numpy as np
from get_epsilon import get_epsilon
from opacity_interpolator import opacity_interpolator
from equation_of_state import equation_of_state
import constants

G = constants.G
nabla_ad = constants.nabla_ad
a_rad = constants.a
c_light = constants.c

def derivs(m, state, X, Y, Z, filename):
    """
    Calculates the derivatives of the four stellar structure equations.
    m: Current enclosed mass (independent variable)
    state: [l, p, r, t] (dependent variables)
    """
    l, p, r, t = state

    # Ensure quantities are not zero

    r = max(r, 1e-5)
    p = max(p, 1e-10)
    t = max(t, 1e-2)

    # 1. Physics
    # Use previously defined EoS, Epsilon, and Opacity functions
    rho = equation_of_state(p, t)
    eps = get_epsilon(rho, t, X, Y, Z)
    kappa = opacity_interpolator(np.log10(rho), np.log10(t), X, Y, Z, filename)
        
    # 2. The First Three ODEs (Mass-based)
    dl_dm = eps
    dp_dm = - (G * m) / (4.0 * np.pi * r**4)
    dr_dm = 1.0 / (4.0 * np.pi * r**2 * rho)

    # 3. Energy Transport (Temperature Gradient)
    
    nabla_rad = (3.0 / (16.0 * np.pi * a_rad * c_light * G)) * (p * kappa * l) / (m * t**4)
    
    # Schwarzschild Criterion: Choose the smaller gradient
    if nabla_rad <= nabla_ad:
        nabla = nabla_rad
    else:
        nabla = nabla_ad
        
    # Temperature ODE
    dt_dm = - (G * m * t) / (4.0 * np.pi * r**4 * p) * nabla

    return np.array([dl_dm, dp_dm, dr_dm, dt_dm])