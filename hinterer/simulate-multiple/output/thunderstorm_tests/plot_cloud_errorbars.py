import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from matplotlib.patches import Ellipse


def compute_para_perp(x_errors, y_errors, azimuths):
    """Compute parallel and perpendicular errors"""
    para_errors = x_errors*np.cos(azimuths) - y_errors*np.sin(azimuths)
    perp_errors = x_errors*np.sin(azimuths) + y_errors*np.cos(azimuths)
    return para_errors, perp_errors


def compute_para_perp_std(x_stds, y_stds, azimuths):
    """Compute standard deviations for parallel and perpendicular components
    using error propagation formula for uncorrelated variables"""
    cos_az = np.cos(azimuths)
    sin_az = np.sin(azimuths)
    
    # Variance of para = (cos(az))^2 * var_x + (sin(az))^2 * var_y
    para_var = (cos_az**2) * (x_stds**2) + (sin_az**2) * (y_stds**2)
    
    # Variance of perp = (sin(az))^2 * var_x + (cos(az))^2 * var_y
    perp_var = (sin_az**2) * (x_stds**2) + (cos_az**2) * (y_stds**2)
    
    return np.sqrt(para_var), np.sqrt(perp_var)


# Get default matplotlib colours
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
dblue = default_colors[0]
dyellow = default_colors[1]
dgreen = default_colors[2]
dred = default_colors[3]
dBlues = LinearSegmentedColormap.from_list('dblue_to_white', [(1, 1, 1), dblue], N=100)
dYellows = LinearSegmentedColormap.from_list('dyellow_to_white', [(1, 1, 1), dyellow], N=100)

# Define file paths for CSV files
file_paths = [
    'fitting_results_hinterer_test_first100.csv',
    'fitting_results_hinterer_test_first100.csv',
    'fitting_results_hinterer.csv',
]

model_names = ['blah1', 'blah2', 'Hinterer on Hinterer']
datasets = []

# Load and process data from CSV files
subdirectory = "./"

for i, file_path in enumerate(file_paths):
    full_path = os.path.join(subdirectory, file_path)
    
    if os.path.exists(full_path):
        # Load data directly from CSV
        df = pd.read_csv(full_path)
        
        # Calculate errors if they're not already in the CSV
        if 'para_err' not in df.columns or 'perp_err' not in df.columns:
            # Convert azimuths to radians for calculation if they're in degrees
            az_rad = df['az_tru'] * np.pi / 180 if df['az_tru'].max() > 6.28 else df['az_tru']
            
            # Calculate parallel and perpendicular errors
            para_errors, perp_errors = compute_para_perp(
                df['x_tru'] - df['x_est'], 
                df['y_tru'] - df['y_est'], 
                az_rad
            )
            
            df['para_err'] = para_errors
            df['perp_err'] = perp_errors
            
            # Calculate standard deviations for para and perp if standard deviations are available
            if 'x_std' in df.columns and 'y_std' in df.columns:
                para_stds, perp_stds = compute_para_perp_std(
                    df['x_std'], 
                    df['y_std'], 
                    az_rad
                )
                
                df['para_std'] = para_stds
                df['perp_std'] = perp_stds
        
        datasets.append(df)
    else:
        print(f"Warning: File not found: {full_path}")
        # Add an empty DataFrame to maintain index alignment
        datasets.append(pd.DataFrame())

for dataset in datasets:
    # Convert rad to deg
    angle_columns = ["inc_tru", "az_tru", "inc_est", "az_est", "inc_err", "az_err"]
    dataset[angle_columns] = dataset[angle_columns] * 180 / np.pi
    
    # Also convert standard deviations if present
    if 'inc_std' in dataset.columns and 'az_std' in dataset.columns:
        std_angle_columns = ["inc_std", "az_std"]
        dataset[std_angle_columns] = dataset[std_angle_columns] * 180 / np.pi

    # Account for wrapping of inc around 180 degrees:
    dataset["inc_err"] = np.mod(dataset["inc_err"], 180)
    dataset["inc_err"] = np.minimum(np.abs(dataset["inc_err"]), 180 - np.abs(dataset["inc_err"]))

    # Account for wrapping of az around 360 degrees:
    dataset["az_err"] = np.mod(dataset["az_err"], 360)
    dataset["az_err"] = np.minimum(np.abs(dataset["az_err"]), 360 - np.abs(dataset["az_err"]))


output_dir = './'

# Function to add an uncertainty ellipse to the plot
def add_uncertainty_ellipse(ax, x, y, x_std, y_std, color='gray', alpha=0.1):
    """
    Add an ellipse to the ax that represents the uncertainty in x and y.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes to add the ellipse to.
    x, y : float
        The center of the ellipse.
    x_std, y_std : float
        The standard deviations in x and y directions.
    color : str, optional
        The color of the ellipse.
    alpha : float, optional
        The transparency of the ellipse.
    """
    # Create 95% confidence ellipse (2 standard deviations)
    ellipse = Ellipse((x, y), width=4*x_std, height=4*y_std, 
                      edgecolor=color, facecolor=color, alpha=alpha)
    ax.add_patch(ellipse)


# Generate plots for each fixed inclination
for inclination in [0, 23, 45, 68, 90]:
   
    # Create a 3x3 figure: 
    # - First column: Theta error colored cloud
    # - Second column: Phi error colored cloud
    # - Third column: Uncertainty cloud
    fig, axs = plt.subplots(3, 3, figsize=(15, 15))

    for i, dataset in enumerate(datasets):

        dataset_inc = dataset[abs(dataset['inc_tru'] - inclination) <= 5]

        if dataset_inc.empty:
          continue

        IQR_multiplier = 100
        
        # --- First column: Error cloud colored by theta error ---
        para_err = dataset_inc["para_err"]
        perp_err = dataset_inc["perp_err"]
        combined_err = np.sqrt(para_err**2 + perp_err**2)
        
        # IQR to remove outliers based on combined error
        Q1 = np.percentile(combined_err, 25)
        Q3 = np.percentile(combined_err, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        
        # Filter data and calculate outlier percentage
        mask = (combined_err >= lower_bound) & (combined_err <= upper_bound)
        total_points = len(combined_err)
        outliers = total_points - np.sum(mask)
        outlier_percent = (outliers / total_points) * 100 if total_points > 0 else 0
        
        # Filter all relevant data using the mask
        filtered_para = para_err[mask]
        filtered_perp = perp_err[mask]
        filtered_inc_err = dataset_inc["inc_err"][mask]  # Using inc_err for coloring
        
        # Single scatter plot with para_err on x-axis and perp_err on y-axis
        scatter = axs[i, 0].scatter(filtered_para, filtered_perp, s=2, c=filtered_inc_err, cmap='rainbow')
 
        # Add axis labels and title
        axs[i, 0].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[i, 0].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[i, 0].set_title(
            f"{model_names[i]}\nθ={inclination}°\nColored by θ error"
        )
        
        # Add colorbar
        cbar = fig.colorbar(scatter, ax=axs[i, 0], pad=0.01)
        cbar.set_label('Theta error (°)')
        
        axs[i, 0].set_xlim(-20, 20)
        axs[i, 0].set_ylim(-20, 20)
        
        # Force aspect ratio to be 1:1
        axs[i, 0].set_aspect('equal')
        axs[i, 0].grid(True, linestyle='--', alpha=0.7)
        axs[i, 0].axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        axs[i, 0].axvline(x=0, color='gray', linestyle='-', alpha=0.5)

        # --- Second column: Error cloud colored by phi error ---
        para_err = dataset_inc["para_err"]
        perp_err = dataset_inc["perp_err"]
        combined_err = np.sqrt(para_err**2 + perp_err**2)
        
        # IQR to remove outliers based on combined error
        Q1 = np.percentile(combined_err, 25)
        Q3 = np.percentile(combined_err, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        
        # Filter data and calculate outlier percentage
        mask = (combined_err >= lower_bound) & (combined_err <= upper_bound)
        
        # Filter all relevant data using the mask
        filtered_para = para_err[mask]
        filtered_perp = perp_err[mask]
        filtered_az_err = dataset_inc["az_err"][mask]  # Using az_err for coloring
        
        # Single scatter plot with para_err on x-axis and perp_err on y-axis
        scatter = axs[i, 1].scatter(filtered_para, filtered_perp, s=2, c=filtered_az_err, cmap='rainbow')
 
        # Add axis labels and title
        axs[i, 1].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[i, 1].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[i, 1].set_title(
            f"{model_names[i]}\nθ={inclination}°\nColored by φ error"
        )
        
        # Add colorbar
        cbar = fig.colorbar(scatter, ax=axs[i, 1], pad=0.01)
        cbar.set_label('Phi error (°)')
        
        axs[i, 1].set_xlim(-20, 20)
        axs[i, 1].set_ylim(-20, 20)
        
        # Force aspect ratio to be 1:1
        axs[i, 1].set_aspect('equal')
        axs[i, 1].grid(True, linestyle='--', alpha=0.7)
        axs[i, 1].axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        axs[i, 1].axvline(x=0, color='gray', linestyle='-', alpha=0.5)

        # --- Third column: Uncertainty visualization ---
        # Only if standard deviations are available
        if 'para_std' in dataset_inc.columns and 'perp_std' in dataset_inc.columns:
            para_err = dataset_inc["para_err"]
            perp_err = dataset_inc["perp_err"]
            para_std = dataset_inc["para_std"]
            perp_std = dataset_inc["perp_std"]
            
            combined_err = np.sqrt(para_err**2 + perp_err**2)
            
            # IQR to remove outliers based on combined error
            Q1 = np.percentile(combined_err, 25)
            Q3 = np.percentile(combined_err, 75)
            IQR = Q3 - Q1
            lower_bound = Q1 - IQR_multiplier * IQR
            upper_bound = Q3 + IQR_multiplier * IQR
            
            # Filter data and calculate outlier percentage
            mask = (combined_err >= lower_bound) & (combined_err <= upper_bound)
            
            # Filter all relevant data using the mask
            filtered_para = para_err[mask]
            filtered_perp = perp_err[mask]
            filtered_para_std = para_std[mask]
            filtered_perp_std = perp_std[mask]
            
            # Compute mean uncertainty
            mean_para_std = np.mean(filtered_para_std)
            mean_perp_std = np.mean(filtered_perp_std)
            
            # Scatter plot but with semi-transparent points
            axs[i, 2].scatter(filtered_para, filtered_perp, s=2, color=dblue, alpha=0.3)
            
            # Sample a subset of points for uncertainty ellipses to avoid clutter
            n_points = len(filtered_para)
            if n_points > 50000000:
                sample_indices = np.random.choice(n_points, 50, replace=False)
            else:
                sample_indices = range(n_points)
                
            # Add uncertainty ellipses for each sampled point
            for j in sample_indices:
                add_uncertainty_ellipse(
                    axs[i, 2], 
                    filtered_para.iloc[j], 
                    filtered_perp.iloc[j],
                    filtered_para_std.iloc[j], 
                    filtered_perp_std.iloc[j],
                    color=dblue,
                    alpha=0.1
                )
            
            # Add a single ellipse showing the average uncertainty
            add_uncertainty_ellipse(
                axs[i, 2], 
                0, 0,  # Center at origin
                mean_para_std, 
                mean_perp_std,
                color=dred,
                alpha=0.2
            )
            
            # Add axis labels and title
            axs[i, 2].set_xlabel('Parallel error (Δ$\parallel$), nm')
            axs[i, 2].set_ylabel('Perpendicular error (Δ$\perp$), nm')
            axs[i, 2].set_title(
                f"{model_names[i]}, θ={inclination}°\nellipses = 95% confidence\n"
                f"Mean σ_para = {mean_para_std:.4f}, Mean σ_perp = {mean_perp_std:.4f}"
            )
            
            axs[i, 2].set_xlim(-20, 20)
            axs[i, 2].set_ylim(-20, 20)
            
            # Force aspect ratio to be 1:1
            axs[i, 2].set_aspect('equal')
            axs[i, 2].grid(True, linestyle='--', alpha=0.7)
            axs[i, 2].axhline(y=0, color='gray', linestyle='-', alpha=0.5)
            axs[i, 2].axvline(x=0, color='gray', linestyle='-', alpha=0.5)
        else:
            # If standard deviations aren't available
            axs[i, 2].text(0.5, 0.5, 'No standard deviation data available', 
                           horizontalalignment='center', verticalalignment='center',
                           transform=axs[i, 2].transAxes)
            axs[i, 2].set_xlim(-20, 20)
            axs[i, 2].set_ylim(-20, 20)
            axs[i, 2].set_aspect('equal')

    plt.tight_layout()
    output_filename = f"cloud_hinterer_inc{round(inclination)}_with_uncertainty.pdf"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()
