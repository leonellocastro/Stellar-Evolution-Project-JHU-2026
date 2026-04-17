import numpy as np

def get_psi(T6):
    # Regime 1: ppI Dominant
    if T6 < 10:
        # linear ramp
        # Transitioning from 1.0 to ~1.3
        psi = 1.0 + 0.03 * (T6 - 1) 
            
    # Regime 2: ppII Dominant
    elif 10 <= T6 < 25:
        # Transitioning from ~1.3 to its peak near 1.95
        psi = 1.3 + 0.043 * (T6 - 10)
            
    # Regime 3: ppIII Dominant
    else:
        # This prevents psi from staying at 2 indefinitely.
        psi = 1.95 - 0.01 * (T6 - 25)
        if psi < 1.45:
            psi = 1.45 # Floor value for pure ppIII
    return psi