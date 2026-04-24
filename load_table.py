import numpy as np
import re

# Find the T, R, and kappa grids for the found table
def load_table(filename, target_id):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    data_block_start = -1
    X, Y, Z = None, None, None

    # 1. Locate the actual data block for Table 1 (skipping the summary)
    for i, line in enumerate(lines):
        # Look for table
        match = re.search(r'TABLE #\s*' + str(target_id) + r'\b', line)
        if match and i > 200:  # Ensure it is the past the summary section
            data_block_start = i
            # EXTRACT METADATA IMMEDIATELY while it is the right line
            X = float(re.search(r'X=([\d\.]+)', line).group(1))
            Y = float(re.search(r'Y=([\d\.]+)', line).group(1))
            Z = float(re.search(r'Z=([\d\.]+)', line).group(1))
            break
    if data_block_start == -1:
        raise ValueError(f"Table {target_id} not found in file.")

    # 2. Extract log R headers from the table's header line
    header_line = lines[data_block_start + 4]
    log_r = [float(x) for x in header_line.split()[1:]]
    
    # 3. Read the opacity rows into lists
    log_T = []
    kappa_matrix = []
    
    for line in lines[data_block_start + 6:]:
        parts = line.split()
        if not parts: continue
        
        try:
            t_val = float(parts[0])
            if not (3.0 <= t_val <= 10.0): break
            
            row = []
            for val_str in parts[1:]:
                v = float(val_str)
                # Convert the 9.999 null to NaN
                row.append(v if v != 9.999 else np.nan)
                
            if len(row) == len(log_r):
                log_T.append(t_val)
                kappa_matrix.append(row)
        except ValueError:
            break
        
    # Return final NumPy arrays
    return np.array(log_T), np.array(log_r), np.array(kappa_matrix), X, Y, Z