import numpy as np
import os, importlib
import test
import time
from save_results import save_results_to_py

def load_existing_results(filename):
    """
    Load existing theta (inc_tru) and phi (az_tru) values from the results file.
    Returns a set of (theta, phi) tuples.
    """
    # Check if the file exists
    if not os.path.exists(filename):
        return set()
    
    # Dynamically load the module
    spec = importlib.util.spec_from_file_location("existing_results", filename)
    existing_results = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(existing_results)
    
    # Extract inc_tru and az_tru lists
    if hasattr(existing_results, 'inc_tru') and hasattr(existing_results, 'az_tru'):
        # Create a set of (theta, phi) tuples
        return {(inc, az) for inc, az in zip(existing_results.inc_tru, existing_results.az_tru)}
    
    return set()

def run_test(phi, theta, filename, num_repeats=20, start_counter=1):
    """
    Runs the Mortensen fit from test.py directly via function call.
    """
    counter = start_counter
    for _ in range(num_repeats):
        try:
            print('-------------------')
            print(f'ϕ = {round(phi*180/np.pi)}°/360°')
            start_time = time.time()
            results, ground_truth = test.run_mortensen_fit(phi, theta)
            results = list(map(float, results))  
            ground_truth = list(map(float, ground_truth))  
#            print(f"Mortensen fitting results are: {results}")
#            print(f"Extracted ground truth string: {ground_truth}")
            # Save results
            save_results_to_py(results, ground_truth, filename)
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f'    {elapsed_time:.2f} seconds')
            
            counter += 1
        except ValueError as e:
            print(f"Error parsing results: {e}")
            continue
    
    return counter

def main():
    # Define theta and phi values
    chosen_thetas = [67.5]#[0, 22.5, 45, 67.5, 90]
    #chosen_phis = [0, 45, 135, 180, 225, 270, 315]
    fixed_phis = np.arange(0, 360, 45)
    #fixed_thetas = list(range(0, 91, 2))
    #done = np.append(np.arange(0, 90, 1), np.arange(90, 361, 5))
    #fixed_phis = [x for x in np.arange(0, 360, 0.5) if x not in done]
    num_repeats = 1
    
    # Initialize counter outside the loop
    counter = 1
    
    filename = f"fitting_results_mortensen_allinone.py"
   
    existing_combinations = load_existing_results(filename)

    # Run tests for each theta and varying phi values
    for theta in chosen_thetas:
        for phi in fixed_phis:

            # Check if this combination already in results file
            if (theta*np.pi/180, phi*np.pi/180) in existing_combinations:
#                print(f"Skipping existing combination: θ = {theta}°, ϕ = {phi}°")
                continue
            
            counter = run_test(phi * np.pi / 180, theta * np.pi / 180, filename, num_repeats, counter)

#    for phi in chosen_phis:
#        filename = f"results_phi_{phi}.csv"
#        for theta in fixed_thetas:
#            counter = run_test(phi * np.pi / 180, theta * np.pi / 180, filename, num_repeats, counter)

if __name__ == "__main__":
    main()
