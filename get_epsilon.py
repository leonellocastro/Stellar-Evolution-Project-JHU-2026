import numpy as np
import constants
from get_f11 import get_f11
from get_psi import get_psi
from get_g11 import get_g11
from get_g14_1 import get_g14_1
from energy_generation import pp_chain, cno_cycle, total_energy_generation

def get_epsilon(rho, T, X, Y, Z):
    T6 = T/1e6
    T7 = T/1e7
    T9 = T/1e9
    components = [X, Y, Z]
    atomic_numbers = [1, 2, 6]      # atomic numbers (H, He, C)
    atomic_masses = [1, 4, 12]      # atomic masses (H, He, C)
    zeta = constants.zeta

    for i in range(len(components)):
        zeta = zeta + atomic_numbers[i]*(atomic_numbers[i] + 1)*components[i]/atomic_masses[i]

    Z1 = constants.Z1      # proton-proton chain
    Z2 = constants.Z2      # proton-proton chain
    
    f11 = get_f11(Z1, Z2, zeta, rho, T7)
    psi = get_psi(T6)
    g11 = get_g11(T9)
    epsilon_pp = pp_chain(psi, f11, g11, rho, X, T9)

    g14_1 = get_g14_1(T9)
    X_CNO = Z
    epsilon_CNO = cno_cycle(g14_1, X_CNO, X, rho, T9)
    epsilon_total = total_energy_generation(epsilon_pp, epsilon_CNO)
    return epsilon_total