import numpy as np
import pandas as pd # Optional, but makes CSV saving easy
import pickle
import constants
from equation_of_state import equation_of_state
from get_epsilon import get_epsilon
from opacity_interpolator import opacity_interpolator
from derivs import derivs
from scipy.integrate import solve_ivp
from shootf import load1, load2

# Load the saved solution and kappa_avg from the previous run
with open('model_result.pkl', 'rb') as f:
    data = pickle.load(f)

sol = data['sol']
kappa_avg = data['kappa_avg']

delta_m = constants.delta_m
filename = constants.filename
M_target = constants.M_target
X = constants.X
Y = constants.Y
Z = constants.Z
a_rad = constants.a
c_light = constants.c
nabla_ad = constants.nabla_ad
G = constants.G

M_mid = M_target / 2.0

Pc, Tc, L_star, R_star = sol.x

y_start_core = load1(delta_m, Pc, Tc)
y_start_surf = load2(M_target, L_star, R_star, kappa_avg)

if sol.success:
    # 1. Re-run integrations with converged values to get the full profile
    # Use dense_output=True to sample at specific mass points
    sol_out = solve_ivp(derivs, [delta_m, M_mid], y_start_core, args=(X, Y, Z, filename), dense_output=True)
    sol_in = solve_ivp(derivs, [M_target, M_mid], y_start_surf, args=(X, Y, Z, filename), dense_output=True)

    # 2. Create a unified mass grid (e.g., 200 points)
    m_points = np.linspace(delta_m, M_target, 200)
    data_rows = []

    for m in m_points:
        # Get [L, P, r, T] from the appropriate solver
        if m <= M_mid:
            state = sol_out.sol(m)
        else:
            state = sol_in.sol(m)
        
        L, P, r, T = state
        rho = equation_of_state(P, T)
        eps = get_epsilon(rho, T, X, Y, Z)
        kappa = opacity_interpolator(np.log10(rho), np.log10(T), X, Y, Z, filename)
        
        # Calculate Nabla (same logic as in derivs)
        nabla_rad = (3.0 / (16.0 * np.pi * a_rad * c_light * G)) * (P * kappa * L) / (m * T**4)
        nabla = min(nabla_rad, nabla_ad)
        nature = "Convective" if nabla_rad > constants.nabla_ad else "Radiative"

        data_rows.append([m, r, rho, T, P, L, eps, kappa, constants.nabla_ad, nabla_rad, nabla, nature])

    # 3. Save to CSV
    cols = ['m', 'r', 'rho', 'T', 'P', 'L', 'epsilon', 'kappa', 'nabla_ad', 'nabla_rad', 'nabla', 'nature']
    df = pd.DataFrame(data_rows, columns=cols)

    # 4. Format to 3 decimal points in scientific notation
    # Apply the format to everything except the 'nature' column
    for col in cols[:-1]:
        df[col] = df[col].map(lambda x: f"{x:.3e}")

    df.to_csv("stellar_structure_5Msun.csv", index=False)
    print("\nMachine-readable table saved to stellar_structure_5Msun.csv")

save_data = {'R_star': R_star, 'L_star': L_star}
with open('star_properties.pkl', 'wb') as f:
    pickle.dump(save_data, f)
print("Star properties saved to star_properties.pkl")