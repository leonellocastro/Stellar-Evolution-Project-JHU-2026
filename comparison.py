import pandas as pd
import numpy as np
import pickle
import constants
from percentage_error import percentage_error
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import cumulative_trapezoid

# 1. Load model
my_model = pd.read_csv('stellar_structure_5Msun.csv')

# 2. Load MESA profile
# Skip the first 5 rows of header, and use whitespace as delimiter
mesa_profile = pd.read_csv('profile.data', skiprows=5, sep=r'\s+')

# 3. Define Constants
M_sun = 1.9884e33
R_sun = 6.957e10
L_sun = 3.828e33

G = constants.G
SIGMA = constants.SIGMA

M_star = constants.M_target

# 4. Prepare MESA data
mesa_profile = mesa_profile.sort_values('mass')

m_mesa_g = mesa_profile['mass'] * M_sun
T_mesa = 10**mesa_profile['logT']
P_mesa = 10**mesa_profile['logP']
rho_mesa = 10**mesa_profile['logRho']
radius_mesa = 10**mesa_profile['logR'] * R_sun

# 5. Calculate Luminosity Profile
# Integrate epsilon over mass to get luminosity
epsilon = mesa_profile['pp'] + mesa_profile['cno'] + mesa_profile['tri_alpha']

# cumulative_trapezoid computes the integral of epsilon dm
L_mesa_erg = cumulative_trapezoid(epsilon, m_mesa_g, initial=0)

# 6. Plotting
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
plt.subplots_adjust(hspace=0.3, wspace=0.3)

# --- GLOBAL FONT SETTINGS ---
plt.rcParams.update({
    'font.size': 20,           # Base font size (affects numbers on axes)
    'axes.titlesize': 20,      # Subplot titles
    'axes.labelsize': 30,      # X and Y labels
    'legend.fontsize': 20,     # Legend font
    'xtick.labelsize': 20,     # Individual axis numbers (X)
    'ytick.labelsize': 20      # Individual axis numbers (Y)
})

# Temperature Plot
axs[0, 0].plot(my_model['m'] / M_star, my_model['T'], 'r-', label='Model')
axs[0, 0].plot(m_mesa_g / M_star, T_mesa, 'k--', label='MESA')
axs[0, 0].set_title('Temperature Profile (Log)')
axs[0, 0].set_ylabel('T (K)', fontsize=24)
axs[0, 0].set_yscale('log')

# Pressure Plot
axs[0, 1].plot(my_model['m'] / M_star, my_model['P'], 'b-', label='Model')
axs[0, 1].plot(m_mesa_g / M_star, P_mesa, 'k--', label='MESA')
axs[0, 1].set_title('Pressure Profile (Log)')
axs[0, 1].set_ylabel(r'$\mathrm{P} \ (\mathrm{dyn} \cdot \mathrm{cm}^{-2})$', fontsize=24)
axs[0, 1].set_yscale('log')

# Radius Plot
# axs[1, 0].plot(my_model['m'] / M_star, my_model['r'], 'g-', label='Model')
# axs[1, 0].plot(m_mesa_g / M_star, radius_mesa, 'k--', label='MESA')
# axs[1, 0].set_title('Radius Profile (Log)')
# axs[1, 0].set_ylabel(r'$\mathrm{R} \ (\mathrm{cm})$')
# axs[1, 0].set_yscale('log')

# Density Plot
axs[1, 0].plot(my_model['m'] / M_star, my_model['rho'], 'g-', label='Model')
axs[1, 0].plot(m_mesa_g / M_star, rho_mesa, 'k--', label='MESA')
axs[1, 0].set_title('Density Profile (Log)')
axs[1, 0].set_ylabel(r'$\rho \ (\mathrm{g} \cdot \mathrm{cm}^{-3})$', fontsize=24)
axs[1, 0].set_yscale('log')

# Luminosity Plot
axs[1, 1].plot(my_model['m'] / M_star, my_model['L'], 'm-', label='Model')
axs[1, 1].plot(m_mesa_g / M_star, L_mesa_erg, 'k--', label='MESA')
axs[1, 1].set_title('Luminosity Profile (Log)')
axs[1, 1].set_ylabel(r'$\mathrm{L} \ (\mathrm{erg} \cdot \mathrm{s}^{-1})$', fontsize=24)
axs[1, 1].set_yscale('log')

for ax in axs.flat:
    ax.set_xlabel(r'M/$M_{\ast}$', fontsize=18)

    ax.yaxis.get_offset_text().set_fontsize(20) 
    ax.xaxis.get_offset_text().set_fontsize(20)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # --- THIS IS THE KEY LINE FOR THE NUMBERS ---
    # axis='both' handles X and Y, labelsize sets the number size
    ax.tick_params(axis='both', which='major', labelsize=20)

plt.suptitle(r'5 $M_{\odot}$ ZAMS: Custom Model vs. MESA', fontsize=30)
plt.savefig('comparison_plot.png')
plt.show()

# Load the star radius and luminosity from the star properties
with open('star_properties.pkl', 'rb') as f:
    data = pickle.load(f)

R_star_model = data['R_star']
L_star_model = data['L_star']

# Calculate model effective temperature using the Stefan-Boltzmann law
T_eff_model = (L_star_model / (4 * np.pi * R_star_model**2 * SIGMA))**0.25

# Calculate model surface gravity
g_surface_model = G * M_star / R_star_model**2

print(f"Model Surface Radius: {R_star_model:.4e} cm")
print(f"Model Surface Radius: {R_star_model/R_sun:.4e} R_sun")
print(f"Model Surface Luminosity: {L_star_model:.4e} erg/s")
print(f"Model Surface Luminosity: {L_star_model/L_sun:.4e} L_sun")
print(f"Model Effective Temperature: {T_eff_model:.4e} K")
print(f"Model Surface Gravity: {g_surface_model:.4e} cm/s^2")

# Compare L, Teff, g and L with MESA values at the surface (last point in the MESA profile)

L_surface_mesa = L_mesa_erg[-1]
R_surface_mesa = radius_mesa.iloc[-1]
T_eff_mesa = (L_surface_mesa / (4 * np.pi * R_surface_mesa**2 * SIGMA))**0.25
g_surface_mesa = G * M_star / R_surface_mesa**2

print(f"MESA Surface Radius: {R_surface_mesa:.4e} cm")
print(f"MESA Surface Radius: {R_surface_mesa/R_sun:.4e} R_sun")
print(f"MESA Surface Luminosity: {L_surface_mesa:.4e} erg/s")
print(f"MESA Surface Luminosity: {L_surface_mesa/L_sun:.4e} L_sun")
print(f"MESA Effective Temperature: {T_eff_mesa:.4e} K")
print(f"MESA Surface Gravity: {g_surface_mesa:.4e} cm/s^2")

# Calculate percentage errors for each quantity
R_error = percentage_error(R_star_model, R_surface_mesa)
L_error = percentage_error(L_star_model, L_surface_mesa)
T_eff_error = percentage_error(T_eff_model, T_eff_mesa)
g_error = percentage_error(g_surface_model, g_surface_mesa)

print(f"Radius Percent Error: {R_error:.4e} %")
print(f"Luminosity Percent Error: {L_error:.4e} %")
print(f"Effective Temperature Percent Error: {T_eff_error:.4e} %")
print(f"Surface Gravity Percent Error: {g_error:.4e} %")

# Plot convective vs radiative zones in the model (nabla_rad vs nabla_ad)

plt.figure(figsize=(10, 6))
model_data = pd.read_csv('stellar_structure_5Msun.csv')
energy_transport = np.where(model_data['nabla_rad'] > model_data['nabla_ad'], 'Convective', 'Radiative')
plt.scatter(model_data['m'] / M_star, model_data['nabla_rad'], c=np.where(energy_transport == 'Convective', 'red', 'blue'), label=r'$\nabla_{rad}$', s=10)
plt.axhline(constants.nabla_ad, color='green', linestyle='--', label=r'$\nabla_{ad}$')
plt.xlabel(r'M/$M_{\ast}$', fontsize=18)
plt.ylabel(r'$\nabla$', fontsize=18)
plt.title('Energy Transport: Convective (Red) vs Radiative (Blue)', fontsize=20)

# Create the legend to show color coding for convective vs radiative zones
proxy_red = Line2D([0], [0], marker='o', color='w', label=r'$\nabla_{\text{rad}}$ (Convective)', markerfacecolor='red', markersize=10)
proxy_blue = Line2D([0], [0], marker='o', color='w', label=r'$\nabla_{\text{rad}}$ (Radiative)', markerfacecolor='blue', markersize=10)
line_ad = Line2D([0], [0], color='green', linestyle='--', label=r'$\nabla_{\text{ad}}$')

# Manually define the legend entries
plt.legend(handles=[proxy_red, proxy_blue, line_ad], fontsize=14, loc='best')
plt.grid(True, alpha=0.3)
plt.savefig('energy_transport.png')
plt.show()

# Find the cutoff mass where the model transitions from convective to radiative
convective_mask = energy_transport == 'Convective'
if np.any(convective_mask):
    cutoff_index = np.where(convective_mask)[0][-1]  # Last convective point
    cutoff_mass = model_data['m'].iloc[cutoff_index]
    print(f"Model transitions from convective to radiative at m = {cutoff_mass:.4e} g ({cutoff_mass/M_star:.4f} M_star)")