import numpy as np
import pandas as pd
import ast
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import warnings
warnings.filterwarnings('ignore')

def count_points_in_confidence_region(points, mean, cov_matrix, confidence_level=0.95):
    """
    Count how many points fall within the specified confidence region.
    
    Parameters:
    -----------
    points : ndarray
        Array of points to check, shape (n_points, n_dimensions)
    mean : ndarray
        Mean vector of the distribution, shape (n_dimensions,)
    cov_matrix : ndarray
        Covariance matrix, shape (n_dimensions, n_dimensions)
    confidence_level : float
        Confidence level (e.g., 0.95 for 95% confidence)
        
    Returns:
    --------
    int
        Number of points within the confidence region
    """
    # Compute the inverse of the covariance matrix
    inv_cov = np.linalg.inv(cov_matrix)
    
    # Degrees of freedom for chi-square (number of dimensions)
    n_dims = len(mean)
    
    # Critical value from chi-square distribution
    # For multivariate normal, squared Mahalanobis distances follow chi-square distribution
    critical_value = stats.chi2.ppf(confidence_level, n_dims)
    
    # Initialize counter
    count = 0
    
    # Check each point
    for point in points:
        # Calculate squared Mahalanobis distance
        diff = point - mean
        mahalanobis_sq = np.dot(np.dot(diff, inv_cov), diff)
        
        # If squared distance is less than critical value, point is within confidence region
        if mahalanobis_sq <= critical_value:
            count += 1
    
    return count

def parse_covariance_matrix(cov_str):
    """Parse the covariance matrix string from the data file."""
    try:
        # Remove outer quotes if present
        if cov_str.startswith('"') and cov_str.endswith('"'):
            cov_str = cov_str[1:-1]
        
        # Handle scientific notation in the string
        cleaned_str = cov_str.replace('e+', 'e').replace('e-', 'e-')
        
        # Parse string to get a numpy array
        cov_matrix = np.array(ast.literal_eval(cleaned_str))
        return cov_matrix
    except Exception as e:
        print(f"Error parsing covariance matrix: {e}")
        # Return identity matrix as fallback
        return np.eye(4)

def is_positive_definite(matrix):
    """Check if a matrix is positive definite."""
    try:
        np.linalg.cholesky(matrix)
        return True
    except np.linalg.LinAlgError:
        return False

def regularize_matrix(matrix):
    """Make a matrix positive definite by adding a small value to the diagonal."""
    # First, ensure the matrix is symmetric
    matrix = 0.5 * (matrix + matrix.T)
    
    # Get eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eigh(matrix)
    
    # Find any non-positive eigenvalues
    min_eig = np.min(eigvals)
    if min_eig <= 0:
        # Add a small positive value to make all eigenvalues positive
        delta = abs(min_eig) + 1e-6
        matrix = matrix + np.eye(len(matrix)) * delta
    
    return matrix

def plot_confidence_ellipse(ax, mean, cov, confidence_level, color='blue', label=None):
    """
    Plot a confidence ellipse on the given axes based on a 2D mean and covariance matrix.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes to plot on
    mean : ndarray
        2D mean vector
    cov : ndarray
        2x2 covariance matrix
    confidence_level : float
        Confidence level (e.g., 0.95 for 95%)
    color : str
        Color for the ellipse
    label : str
        Label for the legend
    """
    # Chi-square value for 2D projection
    chi2_val = stats.chi2.ppf(confidence_level, 2)
    
    # Get eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eigh(cov)
    
    # Sort eigenvalues and eigenvectors
    idx = np.argsort(eigvals)[::-1]  # Sort in descending order
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    
    # Calculate ellipse parameters
    angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))
    width = 2 * np.sqrt(chi2_val * eigvals[0])
    height = 2 * np.sqrt(chi2_val * eigvals[1])
    
    # Create and add the ellipse
    ellipse = Ellipse(
        xy=mean, 
        width=width, 
        height=height, 
        angle=angle,
        edgecolor=color,
        fc='none',
        lw=2,
        label=label
    )
    ax.add_patch(ellipse)

def analyze_dipole(row, dipole_idx, verbose=True, generate_plots=True):
    """
    Analyze a single dipole's confidence regions.
    
    Parameters:
    -----------
    row : pandas.Series
        Row from the dataframe containing the dipole data
    dipole_idx : int
        Dipole index
    verbose : bool
        Whether to print details
    generate_plots : bool
        Whether to generate plots
        
    Returns:
    --------
    dict
        Dictionary of results for this dipole
    """
    try:
        # Parse covariance matrix
        cov_matrix = parse_covariance_matrix(row['covariance'])
        
        # Ensure the matrix is valid
        if not is_positive_definite(cov_matrix):
            if verbose:
                print(f"  Regularizing covariance matrix...")
            cov_matrix = regularize_matrix(cov_matrix)
        
        # Get the parameter estimates (mean vector)
        mean = np.array([
            row['x_est'], 
            row['y_est'], 
            row['inc_est'], 
            row['az_est']
        ])
        
        if verbose:
            print(f"  Mean vector: {mean}")
            print(f"  Covariance matrix shape: {cov_matrix.shape}")
        
        # Print sub-matrices for x-y and inc-az
        if verbose:
            print("\n  Covariance submatrices:")
            print("  x-y:")
            print(cov_matrix[0:2, 0:2])
            print("  inc-az:")
            print(cov_matrix[2:4, 2:4])
        
        # Generate random samples from the distribution
        n_samples = 10000
        np.random.seed(42)
        try:
            points = np.random.multivariate_normal(mean, cov_matrix, size=n_samples)
        except np.linalg.LinAlgError as e:
            print(f"  Error generating samples: {e}")
            print("  Using strongly regularized matrix...")
            cov_matrix = regularize_matrix(cov_matrix)
            cov_matrix = cov_matrix + np.eye(4) * 1e-3  # Add extra regularization
            points = np.random.multivariate_normal(mean, cov_matrix, size=n_samples)
        
        # Check confidence regions
        confidence_levels = [0.68, 0.90, 0.95, 0.99]
        results = {'dipole_index': dipole_idx}
        
        for level in confidence_levels:
            count = count_points_in_confidence_region(points, mean, cov_matrix, level)
            percent = count / n_samples * 100
            results[f"{level*100:.0f}%"] = percent
            
            if verbose:
                print(f"  {level*100:.0f}% confidence region: {percent:.2f}% of points")
        
        # Generate visualizations
        if generate_plots:
            fig, axes = plt.subplots(1, 2, figsize=(16, 8))
            
            # Dimension pairs to plot: x-y and inc-az
            dim_pairs = [
                (0, 1, 'x', 'y', axes[0]),
                (2, 3, 'inclination', 'azimuth', axes[1])
            ]
            
            for dim1, dim2, label1, label2, ax in dim_pairs:
                # Extract the 2D projection
                points_2d = points[:, [dim1, dim2]]
                mean_2d = mean[[dim1, dim2]]
                cov_2d = cov_matrix[np.ix_([dim1, dim2], [dim1, dim2])]
                
                # Print eigenvalues for debugging
                if verbose:
                    eigvals = np.linalg.eigvalsh(cov_2d)
                    print(f"\n  {label1}-{label2} eigenvalues: {eigvals}")
                
                # Plot a subset of points
                n_plot = min(1000, len(points_2d))
                idx = np.random.choice(len(points_2d), n_plot, replace=False)
                ax.scatter(points_2d[idx, 0], points_2d[idx, 1], alpha=0.3, s=5, color='gray')
                
                # Plot the mean
                ax.scatter(mean_2d[0], mean_2d[1], color='red', s=100, marker='x')
                
                # Add ellipses for different confidence levels
                colors = ['blue', 'green', 'orange', 'red']
                for level, color in zip(confidence_levels, colors):
                    plot_confidence_ellipse(
                        ax, 
                        mean_2d, 
                        cov_2d, 
                        level, 
                        color=color, 
                        label=f"{level*100:.0f}%"
                    )
                
                # Set labels and title
                ax.set_xlabel(label1)
                ax.set_ylabel(label2)
                ax.set_title(f"{label1} vs {label2}")
                ax.grid(True, alpha=0.3)
                
                # Add true values if available
                if f"{label1}_tru" in row and f"{label2}_tru" in row:
                    true_val1 = row[f"{label1}_tru"]
                    true_val2 = row[f"{label2}_tru"]
                    ax.scatter(true_val1, true_val2, color='black', s=100, marker='*', label='True Value')
                
                ax.legend()
            
            plt.suptitle(f"Confidence Regions for Dipole {dipole_idx}", fontsize=16)
            plt.tight_layout()
            plt.savefig(f"dipole_{dipole_idx}_confidence_regions.png", dpi=300)
            plt.close()
            
        return results
        
    except Exception as e:
        print(f"Error analyzing dipole {dipole_idx}: {e}")
        import traceback
        traceback.print_exc()
        return {'dipole_index': dipole_idx, 'error': str(e)}

def main():
    # Load the data
    try:
        data = pd.read_csv('fitting_results_hinterer_test.csv')
        print(f"Loaded {len(data)} data points.")
    except Exception as e:
        print(f"Error loading data: {e}")
        return
    
    # Convert numeric columns to float
    for col in data.columns:
        if col != 'covariance':
            try:
                data[col] = pd.to_numeric(data[col], errors='coerce')
            except:
                pass
    
    # Process each dipole
    results = []
    for idx, row in data.iterrows():
        dipole_idx = row['dipole_index']
        print(f"\nAnalyzing dipole {dipole_idx}...")
        
        # Analyze this dipole
        dipole_results = analyze_dipole(row, dipole_idx)
        results.append(dipole_results)
    
    # Create summary
    if results:
        results_df = pd.DataFrame(results)
        print("\nSummary Statistics:")
        print(results_df.describe())
        
        # Plot summary results
        plt.figure(figsize=(12, 8))
        confidence_levels = ['68%', '90%', '95%', '99%']
        for level in confidence_levels:
            if level in results_df.columns:
                plt.plot(results_df['dipole_index'], results_df[level], 'o-', label=f"{level} Confidence")
        
        # Add reference lines
        for level in confidence_levels:
            if level in results_df.columns:
                level_val = float(level.strip('%'))
                plt.axhline(y=level_val, linestyle='--', alpha=0.5)
        
        plt.xlabel('Dipole Index')
        plt.ylabel('Percentage of Points in Confidence Region')
        plt.title('Confidence Region Analysis Results')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig('confidence_region_summary.png')
        plt.close()
        
        # Save results to CSV
        results_df.to_csv('confidence_region_results.csv', index=False)
        print("Results saved to confidence_region_results.csv")

if __name__ == "__main__":
    main()
