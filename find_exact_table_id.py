import numpy as np
import re

# Find table ID by matching X, Y, Z
def find_exact_table_id(filename, target_X, target_Y, target_Z):
    """
    Finds Table ID by matching X, Y, and Z formatted to 4 decimal places.
    """
    # Format targets to match the '0.0000' style in the text file
    sX = f"{target_X:.4f}"
    sY = f"{target_Y:.4f}"
    sZ = f"{target_Z:.4f}"
    
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Search pattern looks for the exact string sequences
    # Example: "X=0.7000 Y=0.2800 Z=0.0200"
    pattern = rf"X={sX}\s*Y={sY}\s*Z={sZ}"

    for line in lines:
        if re.search(pattern, line):
            match_id = re.search(r'TABLE #\s*(\d+)', line)
            if match_id:
                return int(match_id.group(1))
    raise ValueError(f"No table found matching X={sX}, Y={sY}, Z={sZ}")