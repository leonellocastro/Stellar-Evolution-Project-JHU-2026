# shootf.py (load1 and load2 functions)

import numpy as np
from scipy.integrate import solve_ivp
from get_epsilon import get_epsilon
from opacity_interpolator import opacity_interpolator
from equation_of_state import equation_of_state
from load_table import load_table
from derivs import derivs
import constants

# Physical Constants (CGS)
G = constants.G
SIGMA = constants.SIGMA  # Stefan-Boltzmann constant
X = constants.X            # Hydrogen mass fraction
Y = constants.Y         # Helium mass fraction
Z = constants.Z       # Metals mass fraction
m_H = constants.m_H
nabla_ad = constants.nabla_ad
a_rad = constants.a
c_light = constants.c
filename = constants.filename

def load1(delta_m, Pc, Tc):
    """
    Initial values for OUTWARD integration using Mass-based Taylor expansions.
    """
    # 1. Get central physics
    rho_c = equation_of_state(Pc, Tc)
    eps_c = get_epsilon(rho_c, Tc, X, Y, Z)
    kappa_c = opacity_interpolator(np.log10(rho_c), np.log10(Tc), X, Y, Z, filename)
    
    # 2. r and l
    # r = (3m / 4pi*rho)^(1/3)
    r = (3.0 * delta_m / (4.0 * np.pi * rho_c))**(1/3.0)
    # l = epsilon * m
    l = eps_c * delta_m
    
    # 3. Pressure expansion (expressed in M_r)
    # P = Pc - G * (pi/6)^(1/3) * rho^(4/3) * m^(2/3)
    p = Pc - G * (np.pi/6.0)**(1/3.0) * (rho_c**(4/3.0)) * (delta_m**(2/3.0))
    
    # 4. Temperature expansions (expressed in M_r)
    # Calculate central nabla_rad to check for convection
    nabla_rad_c = (3.0 / (16.0 * np.pi * a_rad * c_light * G)) * (Pc * kappa_c * eps_c) / (Tc**4)
    
    if nabla_rad_c <= nabla_ad:
        # Radiative Expansion
        # T^4 = Tc^4 - (1/2ac)*(3/4pi)^(2/3) * kappa * eps * rho^(4/3) * m^(2/3)
        t4 = Tc**4 - (1.0 / (2.0 * a_rad * c_light)) * (3.0/(4.0*np.pi))**(2/3.0) * kappa_c * eps_c * (rho_c**(4/3.0)) * (delta_m**(2/3.0))
        t = t4**0.25
    else:
        # Convective Expansion
        # ln(T) = ln(Tc) - (pi/6)^(1/3) * (G * del_ad * rho^(4/3) / Pc) * m^(2/3)
        ln_t = np.log(Tc) - (np.pi/6.0)**(1/3.0) * (G * nabla_ad * (rho_c**(4/3.0)) / Pc) * (delta_m**(2/3.0))
        t = np.exp(ln_t)
    
    return np.array([l, p, r, t])

def load2(M_star, L_star, R_star, kappa_avg):
    """
    Surface Boundary: Uses average opacity from the composition-specific table.
    """
    T_eff = (L_star / (4.0 * np.pi * R_star**2 * SIGMA))**0.25
    g_surf = (G * M_star) / (R_star**2)
    
    # 3. Surface Pressure (Eddington Boundary)
    P_surf = (2.0 / 3.0) * (g_surf / kappa_avg)
    
    return np.array([L_star, P_surf, R_star, T_eff])

def shootf(guesses, M_star, delta_m, X, Y, Z, kappa_avg, filename):
    if not hasattr(shootf, "iteration"):
        shootf.iteration = 1
    Pc, Tc, L_star, R_star = guesses
    M_mid = M_star / 2.0

    if shootf.iteration > 47:
        print("\n!!! REACHED MAXIMUM ITERATIONS (40) !!!")
        print(f"Last Guesses: Pc={guesses[0]:.2e}, Tc={guesses[1]:.2e}, L={guesses[2]:.2e}, R={guesses[3]:.2e}")
        # Raising an error stops the solver and prevents further computation
        raise StopIteration("Stopping solver: exceeded 47 iterations.")
    
    # 1. Outward Integration (Core to Midpoint)
    y_start_core = load1(delta_m, Pc, Tc)
    sol_out = solve_ivp(derivs, [delta_m, M_mid], y_start_core, args=(X, Y, Z, filename))
    y_mid_out = sol_out.y[:, -1] # Values at M_mid [l, p, r, t]
    
    # 2. Inward Integration (Surface to Midpoint)
    y_start_surf = load2(M_star, L_star, R_star, kappa_avg)
    sol_in = solve_ivp(derivs, [M_star, M_mid], y_start_surf, args=(X, Y, Z, filename))
    y_mid_in = sol_in.y[:, -1] # Values at M_mid [l, p, r, t]
    
    # 3. Calculate the Residual (The "Miss")
    # We want these four differences to eventually be zero

    # 4. Calculate Residuals
    diffs = y_mid_out - y_mid_in
    
    print(f"\nIteration {shootf.iteration}:")

    # Print the "Miss" at the midpoint
    print(f"  > Midpoint L-diff: {diffs[0]:.2e}")
    print(f"  > Midpoint P-diff: {diffs[1]:.2e}")
    print(f"  > Midpoint r-diff: {diffs[2]:.2e}")
    print(f"  > Midpoint T-diff: {diffs[3]:.2e}")

    shootf.iteration += 1

    return diffs