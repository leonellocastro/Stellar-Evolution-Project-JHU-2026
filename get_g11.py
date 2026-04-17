import numpy as np

def get_g11(T9):
    g11 = 1 + 3.82*T9 + 1.51*T9**2 + 0.144*T9**3 - 0.0114*T9**4
    return g11