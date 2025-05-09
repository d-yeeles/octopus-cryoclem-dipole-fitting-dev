import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

def mahalanobis_distance_squared(p_true, p_est, cov):
    """
    Calculate squared Mahalanobis distance between true parameters and estimates.
    
    Parameters:
    - p_true: Ground truth parameters [x, y, theta, phi]
    - p_est: Estimated parameters [x, y, theta, phi] 
    - cov: 4x4 covariance matrix of the parameter estimates
    
    Returns:
    - Squared Mahalanobis distance (scalar)
    """
    try:
        # Make a copy to avoid modifying the originals
        p_true = np.array(p_true, dtype=np.float64).copy()
        p_est = np.array(p_est, dtype=np.float64).copy()
        
        # Handle azimuthal periodicity (phi parameter at index 3)
        # Convert to equivalent angle in [0, 2π)
        p_true[3] = p_true[3] % (2 * np.pi)
        p_est[3] = p_est[3] % (2 * np.pi)
        
        # Choose the smaller angular difference considering periodicity
        phi_diff = p_true[3] - p_est[3]
        if abs(phi_diff) > np.pi:
            if phi_diff > 0:
                phi_diff -= 2 * np.pi
            else:
                phi_diff += 2 * np.pi
        
        # Create the difference vector with the corrected phi difference
        diff = p_true - p_est
        diff[3] = phi_diff
        
        # Handle potential numerical issues in covariance matrix
        # Convert to float64 for better precision
        cov = np.array(cov, dtype=np.float64)
        
        # Ensure the matrix is symmetric
        cov = (cov + cov.T) / 2
        
        # Add a small regularization term to ensure numerical stability
        reg_cov = cov + np.eye(4) * max(1e-8, 1e-10 * np.trace(cov))
        
        # Try standard inversion first with regularized matrix
        try:
            cov_inv = np.linalg.inv(reg_cov)
        except np.linalg.LinAlgError:
            # If that fails, use pseudoinverse
            try:
                cov_inv = np.linalg.pinv(reg_cov)
            except Exception:
                # If all else fails, use diagonal approximation
                cov_inv = np.zeros_like(reg_cov)
                np.fill_diagonal(cov_inv, 1.0 / np.maximum(np.diag(reg_cov), 1e-10))
        
        # Calculate the quadratic form manually to avoid numerical issues
        result = 0.0
        for i in range(4):
            for j in range(4):
                result += diff[i] * cov_inv[i, j] * diff[j]
        
        # Safety check for invalid values
        if not np.isfinite(result) or result < 0:
            # Fallback to diagonal approximation
            diag_result = sum(diff[i]**2 / max(cov[i, i], 1e-10) for i in range(4))
            return diag_result
        
        return result
        
    except Exception as e:
        # Last resort fallback
        print(f"Error in Mahalanobis distance calculation: {e}")
        try:
            # Use a simple normalized Euclidean distance as fallback
            diff = np.array(p_true) - np.array(p_est)
            # Correct azimuthal difference
            if abs(diff[3]) > np.pi:
                diff[3] = min(abs(diff[3]), 2*np.pi - abs(diff[3]))
            # Use diagonal elements as normalization
            variances = np.diag(cov)
            variances = np.maximum(variances, 1e-10)  # Avoid division by zero
            return np.sum(diff**2 / variances)
        except Exception:
            # If everything fails, return a large but finite value
            return 1000.0

def is_within_confidence_region(p_true, p_est, cov, confidence=0.95):
    """
    Test if ground truth is within the confidence region.
    
    Parameters:
    - p_true: Ground truth parameters [x, y, theta, phi]
    - p_est: Estimated parameters [x, y, theta, phi]
    - cov: 4x4 covariance matrix of the parameter estimates
    - confidence: Confidence level (default: 0.95 for 95%)
    
    Returns:
    - Boolean: True if ground truth is within confidence region
    """
    # Calculate squared Mahalanobis distance
    d_squared = mahalanobis_distance_squared(p_true, p_est, cov)
    
    # Critical value from chi-square distribution with 4 degrees of freedom
    critical_value = stats.chi2.ppf(confidence, df=4)
    
    # Check if the point is within the confidence region
    return d_squared <= critical_value

def evaluate_confidence_regions(ground_truths, estimates, covariances, confidence=0.95):
    """
    Evaluate if confidence regions are well-calibrated.
    
    Parameters:
    - ground_truths: List of ground truth parameter vectors
    - estimates: List of parameter estimate vectors
    - covariances: List of covariance matrices
    - confidence: Confidence level to test
    
    Returns:
    - Proportion of ground truths within their confidence regions
    """
    if len(ground_truths) != len(estimates) or len(estimates) != len(covariances):
        raise ValueError("The number of ground truths, estimates, and covariances must be equal")
    
    count_within = 0
    mahalanobis_distances = []
    valid_count = 0
    problem_samples = []
    
    # Get critical value once
    critical_value = stats.chi2.ppf(confidence, df=4)
    print(f"Critical chi-square value for {confidence*100}% confidence with 4 DoF: {critical_value:.4f}")
    
    # Process in smaller batches to avoid memory issues
    batch_size = 50  # Process samples in batches of 50
    total_samples = len(ground_truths)
    
    for batch_start in range(0, total_samples, batch_size):
        batch_end = min(batch_start + batch_size, total_samples)
        print(f"Processing samples {batch_start} to {batch_end-1} of {total_samples}...")
        
        batch_mahalanobis = []
        for i in range(batch_start, batch_end):
            try:
                p_true = ground_truths[i]
                p_est = estimates[i]
                cov = covariances[i]
                
                d_squared = mahalanobis_distance_squared(p_true, p_est, cov)
                mahalanobis_distances.append(d_squared)
                batch_mahalanobis.append(d_squared)
                
                if d_squared <= critical_value:
                    count_within += 1
                elif d_squared > 100:  # Very large value, might indicate a problem
                    problem_samples.append((i, d_squared))
                
                valid_count += 1
                
            except Exception as e:
                print(f"Error processing sample {i}: {e}")
        
        # Write batch results to file to avoid memory buildup
        with open(f'mahalanobis_batch_{batch_start}_{batch_end-1}.txt', 'w') as f:
            for i, d_squared in enumerate(batch_mahalanobis):
                f.write(f"{batch_start + i},{d_squared}\n")
    
    # Write summary statistics
    with open('confidence_region_summary.txt', 'w') as f:
        f.write(f"Total samples: {total_samples}\n")
        f.write(f"Valid samples processed: {valid_count}\n")
        f.write(f"Critical value (95% confidence): {critical_value}\n")
        f.write(f"Samples within confidence region: {count_within}\n")
        f.write(f"Proportion within region: {count_within/valid_count:.4f}\n\n")
        
        # Report statistics on Mahalanobis distances
        if mahalanobis_distances:
            md_array = np.array(mahalanobis_distances)
            f.write("Mahalanobis distance statistics:\n")
            f.write(f"  Min: {np.min(md_array):.4f}\n")
            f.write(f"  Max: {np.max(md_array):.4f}\n")
            f.write(f"  Mean: {np.mean(md_array):.4f}\n")
            f.write(f"  Median: {np.median(md_array):.4f}\n")
            f.write(f"  Std Dev: {np.std(md_array):.4f}\n\n")
        
        # Report problematic samples
        if problem_samples:
            f.write("Potentially problematic samples (MD² > 100):\n")
            for idx, d_squared in problem_samples[:20]:  # List at most 20
                f.write(f"  Sample {idx}: MD² = {d_squared:.4f}\n")
            if len(problem_samples) > 20:
                f.write(f"  ...and {len(problem_samples) - 20} more\n")
    
    print(f"Processed {valid_count} valid samples out of {total_samples} total")
    print(f"Found {count_within} samples within the {confidence*100}% confidence region")
    
    # Create histogram of Mahalanobis distances if we have data
    # But truncate extreme values to make the plot readable
    if mahalanobis_distances:
        try:
            distances_array = np.array(mahalanobis_distances)
            # Clip to reasonable range for visualization
            clipped_distances = np.clip(distances_array, 0, min(100, np.percentile(distances_array, 95) * 2))
            
            plt.figure(figsize=(10, 6))
            plt.hist(clipped_distances, bins=30, alpha=0.7)
            plt.axvline(critical_value, color='r', linestyle='--', 
                      label=f'{confidence*100}% Chi-square critical value')
            plt.xlabel('Squared Mahalanobis Distance (clipped)')
            plt.ylabel('Frequency')
            plt.title('Distribution of Squared Mahalanobis Distances')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig('mahalanobis_distances.png')
            plt.close()  # Close to free memory
            print("Created histogram of Mahalanobis distances")
        except Exception as e:
            print(f"Error creating histogram: {e}")
    
    return count_within / valid_count if valid_count > 0 else 0

def visualize_2d_confidence_ellipse(ax, p_est, cov, indices=(0, 1), confidence=0.95, **kwargs):
    """
    Visualize a 2D projection of the confidence ellipse for two chosen parameters.
    
    Parameters:
    - ax: Matplotlib axis
    - p_est: Estimated parameters 
    - cov: Covariance matrix
    - indices: Tuple of indices to choose which parameters to visualize (default: x, y)
    - confidence: Confidence level
    - **kwargs: Additional arguments for the ellipse
    """
    # Extract the relevant part of the covariance matrix
    cov_2d = np.array([[cov[indices[0], indices[0]], cov[indices[0], indices[1]]],
                       [cov[indices[1], indices[0]], cov[indices[1], indices[1]]]])
    
    # Add a small regularization term for numerical stability
    cov_2d = cov_2d + np.eye(2) * 1e-8
    
    # Find eigenvalues and eigenvectors
    try:
        eigenvals, eigenvecs = np.linalg.eigh(cov_2d)
        
        # Ensure eigenvalues are positive
        eigenvals = np.maximum(eigenvals, 1e-10)
        
        # Chi-square value for 2D with given confidence
        chi2_val = stats.chi2.ppf(confidence, df=2)
        
        # Calculate width and height of ellipse
        width, height = 2 * np.sqrt(chi2_val * eigenvals)
        
        # Calculate angle of ellipse
        angle = np.degrees(np.arctan2(eigenvecs[1, 0], eigenvecs[0, 0]))
        
        # Create and add the ellipse
        ellipse = Ellipse(xy=(p_est[indices[0]], p_est[indices[1]]),
                         width=width, height=height,
                         angle=angle, **kwargs)
        
        ax.add_patch(ellipse)
        return ellipse
    except np.linalg.LinAlgError as e:
        print(f"Error creating ellipse: {e}")
        return None

# Example usage
if __name__ == "__main__":
    # Example: Simulate parameter estimation with ground truth
    np.random.seed(42)
    n_samples = 1000
    
    # Set target confidence level
    confidence_level = 0.95
    
    # Critical chi-square value for 4D
    critical_value = stats.chi2.ppf(confidence_level, df=4)
    print(f"Chi-square critical value for {confidence_level*100}% confidence with 4 DoF: {critical_value:.4f}")
    
    # Lists to store results
    ground_truths = []
    estimates = []
    covariances = []
    
    for i in range(n_samples):
        # Generate random ground truth parameters
        p_true = np.array([np.random.uniform(-10, 10),    # x
                           np.random.uniform(-10, 10),    # y
                           np.random.uniform(-np.pi, np.pi),  # theta
                           np.random.uniform(0, 2*np.pi)])    # phi
        
        # Create a random positive-definite covariance matrix
        A = np.random.randn(4, 4)
        cov = A.T @ A
        # Scale down the covariance to reasonable values
        cov = cov * 0.01
        
        # Generate parameter estimate by adding noise according to the covariance
        noise = np.random.multivariate_normal(np.zeros(4), cov)
        p_est = p_true + noise
        
        ground_truths.append(p_true)
        estimates.append(p_est)
        covariances.append(cov)
    
    # Evaluate if confidence regions are well-calibrated
    proportion_within = evaluate_confidence_regions(ground_truths, estimates, covariances, confidence_level)
    print(f"Proportion of ground truths within {confidence_level*100}% confidence regions: {proportion_within:.4f}")
    
    # Visualize a few examples (2D projections)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    param_names = ['x', 'y', 'theta', 'phi']
    
    # Projections to visualize
    projections = [(0, 1), (0, 2), (1, 3), (2, 3)]
    
    for i, (ax, proj) in enumerate(zip(axes.flatten(), projections)):
        idx1, idx2 = proj
        
        # Plot first 10 examples
        for j in range(10):
            # Plot estimate point
            ax.plot(estimates[j][idx1], estimates[j][idx2], 'bo', alpha=0.5)
            
            # Plot ground truth point
            ax.plot(ground_truths[j][idx1], ground_truths[j][idx2], 'ro')
            
            # Plot confidence ellipse
            visualize_2d_confidence_ellipse(ax, estimates[j], covariances[j], 
                                           indices=proj, confidence=confidence_level,
                                           alpha=0.3, facecolor='blue', edgecolor='blue')
        
        ax.set_xlabel(param_names[idx1])
        ax.set_ylabel(param_names[idx2])
        ax.set_title(f"{param_names[idx1]}-{param_names[idx2]} Projection")
        ax.grid(True)
    
    plt.tight_layout()
    plt.savefig('confidence_regions.png')
    plt.show()
