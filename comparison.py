import pandas as pd
import numpy as np
import constants
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

# 1. Load model
my_model = pd.read_csv('stellar_structure_5Msun.csv')

# 2. Load MESA profile
# Skip the first 5 rows of header, and use whitespace as delimiter
mesa_profile = pd.read_csv('profile.data', skiprows=5, sep=r'\s+')

# 3. Define Constants (MESA usually uses these exact CGS values)
M_sun = 1.9884e33
R_sun = 6.957e10
L_sun = 3.828e33

M_star = constants.M_target

# 4. Prepare MESA data
mesa_profile = mesa_profile.sort_values('mass')

m_mesa_g = mesa_profile['mass'] * M_sun
T_mesa = 10**mesa_profile['logT']
P_mesa = 10**mesa_profile['logP']
rho_mesa = 10**mesa_profile['logRho']

# 5. Calculate Luminosity Profile
# Integrate epsilon over mass to get luminosity
epsilon = mesa_profile['pp'] + mesa_profile['cno'] + mesa_profile['tri_alpha']

# cumulative_trapezoid computes the integral of epsilon dm
L_mesa_erg = cumulative_trapezoid(epsilon, m_mesa_g, initial=0)

# 6. Plotting
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
plt.subplots_adjust(hspace=0.3, wspace=0.3)

# Temperature Plot
axs[0, 0].plot(my_model['m'] / M_star, my_model['T'], 'r-', label='Model')
axs[0, 0].plot(m_mesa_g / M_star, T_mesa, 'k--', label='MESA')
axs[0, 0].set_title('Temperature Profile')
axs[0, 0].set_ylabel('T (K)')

# Pressure Plot
axs[0, 1].plot(my_model['m'] / M_star, my_model['P'], 'b-', label='Model')
axs[0, 1].plot(m_mesa_g / M_star, P_mesa, 'k--', label='MESA')
axs[0, 1].set_title('Pressure Profile (Log)')
axs[0, 1].set_ylabel(r'$\mathrm{P} \ (\mathrm{dyn} \cdot \mathrm{cm}^{-2})$')
axs[0, 1].set_yscale('log')

# Density Plot
axs[1, 0].plot(my_model['m'] / M_star, my_model['rho'], 'g-', label='Model')
axs[1, 0].plot(m_mesa_g / M_star, rho_mesa, 'k--', label='MESA')
axs[1, 0].set_title('Density Profile (Log)')
axs[1, 0].set_ylabel(r'$\mathrm{Rho} \ (\mathrm{g} \cdot \mathrm{cm}^{-3})$')
axs[1, 0].set_yscale('log')

# Luminosity Plot
axs[1, 1].plot(my_model['m'] / M_star, my_model['L'], 'm-', label='Model')
axs[1, 1].plot(m_mesa_g / M_star, L_mesa_erg, 'k--', label='MESA')
axs[1, 1].set_title('Luminosity Profile')
axs[1, 1].set_ylabel(r'$\mathrm{L} \ (\mathrm{erg} \cdot \mathrm{s}^{-1})$')

for ax in axs.flat:
    # set x-axis as a ratio of M/M_star for better comparison
    ax.set_xlabel(r'M/$M_{\ast}$')
    # ax.set_xlabel('Mass (g)')
    ax.legend()
    ax.grid(True, alpha=0.3)

plt.suptitle(r'5 $M_{\odot}$ ZAMS: Custom Model vs. MESA', fontsize=16)
plt.savefig('comparison_plot.png')
plt.show()