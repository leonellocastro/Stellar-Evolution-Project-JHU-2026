import numpy as np

# Calculate percentage error between model and MESA profiles
def percentage_error(model, mesa):
    """
    Calculate percentage error between model and MESA profiles.
    
    Parameters:
    model (array-like): The values from the model.
    mesa (array-like): The corresponding values from MESA.
    
    Returns:
    array-like: Percentage error at each point.
    """
    error = (np.abs(model - mesa) / np.abs(mesa)) * 100.0
    
    return error