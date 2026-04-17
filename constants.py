# constants.py

G = 6.674e-8
SIGMA = 5.670e-5  # Stefan-Boltzmann constant
MSUN = 1.989e33
X = 0.70            # Hydrogen mass fraction
Y = 0.28            # Helium mass fraction
Z = 1 - X - Y       # Metals mass fraction
k = 1.381e-16 # Boltzmann constant in CGS
c = 2.998e10 # Speed of light in cm/s
a = 7.5657e-15 # radiation constant in CGS
N_A = 6.022e23 # Avogadro's number
k_B = 1.38e-16
m_H = 1.67e-24
mu = (2*X + (3/4)*Y + (1/2)*Z)**(-1)
M_target = 5.0 * MSUN
delta_m = 1e-6 * M_target # Start 1 millionth of the mass away from center
Z1 = 1      # proton-proton chain
Z2 = 1      # proton-proton chain
zeta = 0.0 # Initial zeta value for get_epsilon
nabla_ad = 0.4
LSUN = 3.828e33
RSUN = 6.957e10
Pc_guess = 6.8e16
Tc_guess = 2.6e7
L_guess = 550*LSUN
R_guess = 2.5*RSUN
filename = r"opacity_table.txt"