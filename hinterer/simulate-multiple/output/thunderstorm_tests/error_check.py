import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import ast
import sys
import traceback
from error_check_funcs import evaluate_confidence_regions, visualize_2d_confidence_ellipse

# Set matplotlib to use a non-interactive backend to prevent GUI-related crashes
import matplotlib
matplotlib.use('Agg')

# Load the data
try:
    data = pd.read_csv('fitting_results_hinterer_test.csv')
    print(f"Loaded {len(data)} samples from fitting_results_hinterer_test.csv")
except Exception as e:
    print(f"Error loading data: {e}")
    sys.exit(1)

# Function to parse covariance matrix from string
def parse_covariance(cov_str):
    # Clean the string and parse it
    try:
        # First try standard JSON parsing
        cov_matrix = json.loads(cov_str)
    except json.JSONDecodeError:
        try:
            # If that fails, use ast.literal_eval which can handle scientific notation
            cov_matrix = ast.literal_eval(cov_str)
        except Exception as e:
            print(f"Error parsing covariance matrix: {e}")
            # Return a default diagonal matrix as fallback
            return np.eye(4) * 0.1
    
    # Ensure the matrix is 4x4
    if len(cov_matrix) != 4 or any(len(row) != 4 for row in cov_matrix):
        print("Warning: Covariance matrix does not have 4x4 shape. Using default.")
        return np.eye(4) * 0.1
    
    # Convert to numpy array
    cov_array = np.array(cov_matrix, dtype=float)
    
    # Check for invalid values
    if not np.all(np.isfinite(cov_array)):
        print("Warning: Covariance matrix contains NaN or infinite values. Using default.")
        return np.eye(4) * 0.1
    
    # Ensure matrix is symmetric
    cov_array = (cov_array + cov_array.T) / 2
    
    # Add a small regularization to ensure numerical stability
    cov_array = cov_array + np.eye(4) * 1e-8
    
    return cov_array

# Prepare the data in the format required by the script
ground_truths = []
estimates = []
covariances = []
sample_indices = []

for i, row in data.iterrows():
    try:
        # Ground truth parameters [x, y, theta, phi]
        # Here theta corresponds to inclination and phi to azimuth
        p_true = [row.x_tru, row.y_tru, row.inc_tru, row.az_tru]
        
        # Estimated parameters [x, y, theta, phi]
        p_est = [row.x_est, row.y_est, row.inc_est, row.az_est]
        
        # Parse and add covariance matrix
        cov_matrix = parse_covariance(row.covariance)
        
        # Add to our lists
        ground_truths.append(p_true)
        estimates.append(p_est)
        covariances.append(cov_matrix)
        sample_indices.append(i)
        
    except Exception as e:
        print(f"Skipping sample {i} due to error: {e}")

print(f"Successfully processed {len(ground_truths)} out of {len(data)} samples")

# Set target confidence level
confidence_level = 0.95

# Evaluate if confidence regions are well-calibrated
try:
    proportion_within = evaluate_confidence_regions(ground_truths, estimates, covariances, confidence_level)
    print(f"Proportion of ground truths within {confidence_level*100}% confidence regions: {proportion_within:.4f}")
except Exception as e:
    print(f"Error during confidence region evaluation: {e}")
    # Continue with visualization

# Visualize examples (2D projections) - with memory management
try:
    print("\nCreating visualization...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    param_names = ['x', 'y', 'inclination', 'azimuth']

    # Projections to visualize
    projections = [(0, 1), (0, 2), (1, 3), (2, 3)]

    for i, (ax, proj) in enumerate(zip(axes.flatten(), projections)):
        idx1, idx2 = proj
        print(f"  Processing projection {param_names[idx1]}-{param_names[idx2]}...")
        
        # Plot available examples (up to 10)
        for j in range(min(10, len(estimates))):
            try:
                # Plot estimate point
                ax.plot(estimates[j][idx1], estimates[j][idx2], 'bo', alpha=0.5)
                
                # Plot ground truth point
                ax.plot(ground_truths[j][idx1], ground_truths[j][idx2], 'ro')
                
                # Plot confidence ellipse - with explicit error handling
                try:
                    ellipse = visualize_2d_confidence_ellipse(ax, estimates[j], covariances[j], 
                                                  indices=proj, confidence=confidence_level,
                                                  alpha=0.3, facecolor='blue', edgecolor='blue')
                    if ellipse is None:
                        print(f"    Could not create ellipse for sample {j}")
                except Exception as e:
                    print(f"    Error creating ellipse for sample {j}: {e}")
            except Exception as e:
                print(f"    Error plotting sample {j} on subplot {i}: {e}")
                continue
        
        ax.set_xlabel(param_names[idx1])
        ax.set_ylabel(param_names[idx2])
        ax.set_title(f"{param_names[idx1]}-{param_names[idx2]} Projection")
        ax.grid(True)

    plt.tight_layout()
    print("  Saving figure...")
    plt.savefig('confidence_regions_from_data.png')
    print("  Figure saved. Closing figure to free memory.")
    plt.close(fig)  # Close the figure to free memory

    print("Analysis complete. Visualization saved as 'confidence_regions_from_data.png'")
except Exception as e:
    print(f"Error during visualization: {e}")

# Add a detailed analysis of the azimuth values
try:
    print("\nPerforming detailed analysis of azimuth values...")
    az_true = [gt[3] for gt in ground_truths]
    az_est = [est[3] for est in estimates]
    
    # First, save numerical statistics to text file without visualization
    with open('azimuth_stats.txt', 'w') as f:
        f.write("Azimuth Analysis Statistics\n")
        f.write("==========================\n\n")
        
        # Basic statistics
        az_diffs = [abs(t - e) for t, e in zip(az_true, az_est)]
        az_diffs_periodic = []
        
        # Convert estimated azimuths to be closer to the true values (considering periodicity)
        az_est_corrected = []
        for true, est in zip(az_true, az_est):
            # Convert to equivalent angle in [0, 2π)
            true_mod = true % (2 * np.pi)
            est_mod = est % (2 * np.pi)
            
            # Calculate difference considering periodicity
            diff = true_mod - est_mod
            if abs(diff) > np.pi:
                if diff > 0:
                    diff = diff - 2 * np.pi
                else:
                    diff = diff + 2 * np.pi
            
            az_diffs_periodic.append(abs(diff))
            
            # Store corrected estimate for plotting
            if diff > 0:
                est_mod += 2 * np.pi
            elif diff < 0:
                est_mod -= 2 * np.pi
            az_est_corrected.append(est_mod)
        
        # Write statistics
        f.write(f"Number of samples: {len(az_true)}\n")
        f.write(f"Average raw azimuth difference: {np.mean(az_diffs):.4f} radians\n")
        f.write(f"Average periodic azimuth difference: {np.mean(az_diffs_periodic):.4f} radians\n")
        f.write(f"Maximum periodic difference: {np.max(az_diffs_periodic):.4f} radians\n")
        f.write(f"Number of samples with diff > π/2: {sum(1 for d in az_diffs_periodic if d > np.pi/2)}\n")
        
        # Write first 20 sample data
        f.write("\nSample data (first 20 samples):\n")
        f.write("Index  True(rad)  Est(rad)  Raw Diff  Periodic Diff\n")
        for i in range(min(20, len(az_true))):
            f.write(f"{sample_indices[i]:5d}  {az_true[i]:9.4f}  {az_est[i]:8.4f}  {az_diffs[i]:8.4f}  {az_diffs_periodic[i]:12.4f}\n")
    
    print("  Azimuth statistics saved to azimuth_stats.txt")
    
    # Now create the visualization in separate try block to isolate potential errors
    try:
        print("  Creating azimuth visualizations...")
        fig = plt.figure(figsize=(12, 10))
        
        # First subplot: Raw values scatterplot
        plt.subplot(2, 1, 1)
        plt.scatter(sample_indices[:100], az_true[:100], c='r', label='True')
        plt.scatter(sample_indices[:100], az_est[:100], c='b', label='Estimated')
        plt.xlabel('Sample Index')
        plt.ylabel('Azimuth (radians)')
        plt.title('Raw Azimuth Values (first 100 samples)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Second subplot: Histogram of differences
        plt.subplot(2, 1, 2)
        plt.hist(az_diffs_periodic, bins=30, alpha=0.7)
        plt.axvline(np.pi/2, color='r', linestyle='--', label='π/2')
        plt.xlabel('Periodic Azimuth Difference (radians)')
        plt.ylabel('Frequency')
        plt.title('Distribution of Azimuth Differences')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('azimuth_analysis.png')
        plt.close(fig)  # Close the figure to free memory
        
        print("  Azimuth visualization saved as 'azimuth_analysis.png'")
    except Exception as e:
        print(f"  Error during azimuth visualization: {e}")
        traceback.print_exc()
        
    print("Azimuth analysis complete.")
except Exception as e:
    print(f"Error during azimuth analysis: {e}")
    traceback.print_exc()
