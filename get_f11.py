import numpy as np

def get_f11(Z1, Z2, zeta, rho, T7):
    f11 = np.exp((5.92*10e-3)*Z1*Z2*(zeta*rho/(T7**3))**(1/2))
    return f11