import numpy as np
from scipy.interpolate import RegularGridInterpolator
from find_exact_table_id import find_exact_table_id
from load_table import load_table

# Write the interpolation function
def opacity_interpolator(log_rho, log_T, X, Y, Z, filename):
    # 1. Locate the correct table
    try:
        table_id = find_exact_table_id(filename, X, Y, Z)
    except ValueError as e:
        return np.nan

    # 2. Load the data using existing load_table function
    lT_grid, lR_grid, kappa_mat, x_out, y_out, z_out = load_table(filename, table_id)

    # 3. Coordinate Conversion
    # log10(R) = log10(rho) - 3*log10(T) + 18 (simplified version of the T6 version converting rho to R)
    log_R_target = log_rho - 3 * log_T + 18

    # 4. Set up the 2D Interpolator
    # method='linear'
    interp = RegularGridInterpolator(
        (lT_grid, lR_grid), 
        kappa_mat, 
        method='linear', 
        bounds_error=False, 
        fill_value=np.nan
    )

    # 5. Execute Interpolation
    # Point must be [[logT, logR]]
    point = np.array([[log_T, log_R_target]])
    result = interp(point)
    kappa = 10**result[0]  # Convert log10(kappa) back to kappa
    return kappa