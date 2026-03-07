import numpy as np
import constants as const
# Equation of state in terms of P, T, X, Y

def rho(P,T,X,Y):
    # P = P_gas + P_rad
    # P = (rho * k * T) / (mu(X,Y) * m_H) + (1/3) * a * T**4
    # Approximate mean molecular weight for a fully ionized gas
    Z = 1 - X - Y
    mu = (2*X + (3/4)*Y + (1/2)*Z)**(-1)
    rho = (P - (1/3) * const.a * T**4) * mu / (const.N_A * const.k * T)
    return rho

# Examples

X = np.array([0, 0.70])
Y = np.array([0.98, 0.28])
T = np.array([10**7.55, 10**6.91]) # K
P = np.array([10**16.85, 10**16.87]) # dyn/cm^2
mu = (2*X + (3/4)*Y + (1/2)*(1-X-Y))**(-1)

for i in range(len(X)):
    print(f"X = {X[i]:.2f}, Y = {Y[i]:.2f}, T = {T[i]:.2e} K, P = {P[i]:.2e} dyn/cm^2, rho = {rho(P[i], T[i], X[i], Y[i]):.2e} g/cm^3")

rho_value = np.array([rho(P[i], T[i], X[i], Y[i]) for i in range(len(X))])

# Calculate beta = P_gas / P_total

P_gas = (rho_value * const.N_A * const.k * T) / mu

beta = P_gas / P

print(mu)
print(P_gas)
print(beta)