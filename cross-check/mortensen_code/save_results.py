import csv
import os
import ast

def save_results_to_py(results, ground_truth, filename="fitting_results.py"):
    """
    Collects results and ground truth values, and saves them to a Python file
    with lists for true values, estimated values, and errors.
    
    Parameters:
    - results: Tuple of estimated values (phi, theta, x, y, photons)
    - ground_truth: Tuple of ground truth values (phi, theta, x, y, photons)
    - filename: Python filename (default: "fitting_results.py")
    """
    # Extract values from tuples
    phi_est, theta_est, x_est, y_est, photons_est = results
    phi_true, theta_true, x_true, y_true, photons_true = ground_truth
    
    # Calculate errors
    x_err = x_est - x_true
    y_err = y_est - y_true
    inc_err = theta_est - theta_true
    az_err = phi_est - phi_true
    photon_err = photons_est - photons_true
    
    # Define all required keys
    required_keys = [
        'x_tru', 'y_tru', 'inc_tru', 'az_tru',
        'x_est', 'y_est', 'inc_est', 'az_est',
        'x_err', 'y_err', 'inc_err', 'az_err',
        'photon_tru', 'photon_est', 'photon_err'
    ]
    
    # Initialize data dictionary with empty lists for all keys
    data = {key: [] for key in required_keys}
    
    # Read existing data if file exists
    if os.path.isfile(filename):
        try:
            with open(filename, 'r') as file:
                content = file.readlines()
            
            # Parse existing lists
            for line in content:
                if '=' in line:
                    var_name, var_value = line.strip().split('=', 1)
                    var_name = var_name.strip()
                    try:
                        # Safely evaluate the list
                        data[var_name] = ast.literal_eval(var_value.strip())
                    except (SyntaxError, ValueError, KeyError):
                        # Keep as empty list if parsing fails
                        pass
        except Exception as e:
            print(f"Warning: Error reading existing file: {e}")
            # Continue with empty lists if there's any issue
    
    # Ensure all required keys exist
    for key in required_keys:
        if key not in data:
            data[key] = []
    
    # Append new values to lists
    data['x_tru'].append(x_true)
    data['y_tru'].append(y_true)
    data['inc_tru'].append(theta_true)
    data['az_tru'].append(phi_true)
    data['x_est'].append(x_est)
    data['y_est'].append(y_est)
    data['inc_est'].append(theta_est)
    data['az_est'].append(phi_est)
    data['x_err'].append(x_err)
    data['y_err'].append(y_err)
    data['inc_err'].append(inc_err)
    data['az_err'].append(az_err)
    data['photon_tru'].append(photons_true)
    data['photon_est'].append(photons_est)
    data['photon_err'].append(photon_err)
    
    # Write updated lists back to file
    with open(filename, 'w') as file:
        for var_name in required_keys:
            file.write(f"{var_name} = {data[var_name]}\n")
