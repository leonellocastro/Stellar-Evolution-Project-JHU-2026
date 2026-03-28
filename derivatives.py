import numpy as np
import constants as const
import equation_of_state
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

def derivs(m, y, X, Y):
    """
    y[0]=r, y[1]=P, y[2]=L, y[3]=T
    
    """
    r, P, l, T = y
    
    # 1. Get local physics
    rho = equation_of_state(P, T, X, Y) # [cite: 6]
    kappa = get_opacity(rho, T, X, Y) # [cite: 11]
    eps = get_energy_gen(rho, T, X, Y) # [cite: 13]
    
    # 2. Calculate gradients
    # Need Nabla_rad and Nabla_ad to determine if convective 
    # Nabla = min(Nabla_rad, Nabla_ad) for simplicity on ZAMS
    
    # 3. Define the derivatives [cite: 27, 28]
    dr_dm = 1.0 / (4.0 * np.pi * r**2 * rho)
    dP_dm = -(const.G * m) / (4.0 * np.pi * r**4)
    dl_dm = eps

    Nabla = const.Nabla_ad
    Nabla_rad = (3/(16*np.pi*const.a*const.c))*(P*kappa/(T^4))*l/(const.G*m)
    if (Nabla_rad <= const.Nabla_ad):
        Nabla = Nabla_rad
    else:
        Nabla = const.Nabla_ad

    dT_dm = -((const.G*m*T)/(4.0*np.pi*(r**4)* P))*Nabla
    
    return [dr_dm, dP_dm, dl_dm, dT_dm]