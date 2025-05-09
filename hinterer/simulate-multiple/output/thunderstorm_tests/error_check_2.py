import pandas as pd
import numpy as np
import ast
from scipy import stats

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





# Read the CSV file
df = pd.read_csv('fitting_results_hinterer_test.csv')

# First, let's check what columns are actually in the dataframe
print("Columns in the dataframe:", df.columns.tolist())

# Assuming one of these columns contains your matrix data
# Replace 'actual_matrix_column_name' with the correct column name from the printed list
matrix_column_name = 'actual_matrix_column_name'  # Update this after seeing the columns

# Check if the column exists before proceeding
if matrix_column_name in df.columns:
    # Convert the string representation to actual arrays
    df['matrix_numpy'] = df[matrix_column_name].apply(lambda x: np.array(ast.literal_eval(x)))
    
    # Now you can access the matrix as a numpy array
    # For example, to get the first row's matrix:
    matrix = df.loc[0, 'matrix_numpy']
    print("Matrix from first row:")
    print(matrix)
    print("Type:", type(matrix))
else:
    print(f"Column '{matrix_column_name}' not found in the dataframe.")
    # Optionally, show a sample of the dataframe to help identify the right column
    print("\nFirst few rows of the dataframe:")
    print(df.head())






# Mean vector (estimated parameters)
mean = np.array([10.0, 5.0, 45.0, 30.0])  # x, y, inclination, azimuth
    
# Generate some random test points
np.random.seed(42)  # For reproducibility
n_points = 1000
points = np.random.multivariate_normal(mean, cov_matrix, size=n_points)
    
# Compute how many points are within the 95% confidence region
count = count_points_in_confidence_region(points, mean, cov_matrix, confidence_level=0.95)
    
print(f"Number of points within 95% confidence region: {count}")
print(f"Percentage: {count/n_points*100:.2f}%")
