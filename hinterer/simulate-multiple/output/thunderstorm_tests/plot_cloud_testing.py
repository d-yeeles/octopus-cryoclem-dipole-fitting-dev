import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns


def compute_para_perp(x_diffs, y_diffs, azimuths):
    """Compute parallel and perpendicular errors"""
    para_diffs = x_diffs*np.cos(azimuths) - y_diffs*np.sin(azimuths)
    perp_diffs = x_diffs*np.sin(azimuths) + y_diffs*np.cos(azimuths)
    return para_diffs, perp_diffs


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
    './fitting_results_hinterer_test.csv',
    './fitting_results_hinterer_test.csv',
]

model_names = ['hinterer', 'mortensen']
datasets = []

# Load and process data from CSV files
subdirectory = "./"

for i, file_path in enumerate(file_paths):
    full_path = os.path.join(subdirectory, file_path)
    
    if os.path.exists(full_path):
        # Load data directly from CSV
        df = pd.read_csv(full_path)
        
        # Calculate errors if they're not already in the CSV
        if 'para_dif' not in df.columns or 'perp_dif' not in df.columns:
            # Convert azimuths to radians for calculation if they're in degrees
            az_rad = df['az_thu'] * np.pi / 180 if df['az_thu'].max() > 6.28 else df['az_thu']
            
            # Calculate parallel and perpendicular errors
            para_diffs, perp_diffs = compute_para_perp(
                df['x_thu'] - df['x_est'], 
                df['y_thu'] - df['y_est'], 
                az_rad
            )
            
            df['para_dif'] = para_diffs
            df['perp_dif'] = perp_diffs
        
        datasets.append(df)
    else:
        print(f"Warning: File not found: {full_path}")
        # Add an empty DataFrame to maintain index alignment
        datasets.append(pd.DataFrame())

# Process datasets
for dataset in datasets:
    if not dataset.empty:
        # Check if angles are in radians and convert to degrees if needed
        for col in ['inc_thu', 'az_thu', 'inc_est', 'az_est', 'inc_dif', 'az_dif']:
            if col in dataset.columns:
                # If max value is less than 7, assume radians and convert to degrees
                if dataset[col].abs().max() < 7:
                    dataset[col] = dataset[col] * 180 / np.pi

        # Account for wrapping of inclination around 180 degrees
        if 'inc_dif' in dataset.columns:
            dataset['inc_dif'] = np.mod(dataset['inc_dif'], 180)
            dataset['inc_dif'] = np.minimum(np.abs(dataset['inc_dif']), 180 - np.abs(dataset['inc_dif']))

        # Account for wrapping of azimuth around 360 degrees
        if 'az_dif' in dataset.columns:
            dataset['az_dif'] = np.mod(dataset['az_dif'], 360)
            dataset['az_dif'] = np.minimum(np.abs(dataset['az_dif']), 360 - np.abs(dataset['az_dif']))

# Output directory
output_dir = './'

# Generate plots for each fixed inclination
for inclination in [0]:
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    for i, dataset in enumerate(datasets):
        if dataset.empty:
            continue

        # Filter data for given inclination
        dataset_inc = dataset[abs(dataset['inc_thu'] - inclination) <= 5]

        if dataset_inc.empty:
            continue

        IQR_multiplier = 500000000000  # Effectively no outlier removal
        
        # Get error data
        para_dif = dataset_inc["para_dif"]
        perp_dif = dataset_inc["perp_dif"]
        combined_dif = np.sqrt(para_dif**2 + perp_dif**2)
        
        # IQR to remove outliers based on combined error
        Q1 = np.percentile(combined_dif, 25)
        Q3 = np.percentile(combined_dif, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        
        # Filter data and calculate outlier percentage
        mask = (combined_dif >= lower_bound) & (combined_dif <= upper_bound)
        total_points = len(combined_dif)
        outliers = total_points - np.sum(mask)
        outlier_percent = (outliers / total_points) * 100 if total_points > 0 else 0
        
        # Plot inclination error (theta error)
        filtered_para = para_dif[mask]
        filtered_perp = perp_dif[mask]
        filtered_inc_dif = dataset_inc["inc_dif"][mask]
        
        scatter = axs[i, 0].scatter(filtered_para, filtered_perp, s=2, c=filtered_inc_dif, cmap='rainbow')
        
        axs[i, 0].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[i, 0].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[i, 0].set_title(f"{model_names[i]}\nθ={inclination}°")
        
        cbar = fig.colorbar(scatter, ax=axs[i, 0], pad=0.01)
        cbar.set_label('Theta error, °')
        
#        axs[i, 0].set_xlim(-20, 20)
#        axs[i, 0].set_ylim(-20, 20)
        axs[i, 0].set_aspect('equal')
        axs[i, 0].grid(True, linestyle='--', alpha=0.7)
        axs[i, 0].axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        axs[i, 0].axvline(x=0, color='gray', linestyle='-', alpha=0.5)

        # Plot azimuth error (phi error)
        filtered_az_dif = dataset_inc["az_dif"][mask]
        
        scatter = axs[i, 1].scatter(filtered_para, filtered_perp, s=2, c=filtered_az_dif, cmap='rainbow')
        
        axs[i, 1].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[i, 1].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[i, 1].set_title(f"{model_names[i]}\nθ={inclination}°")
        
        cbar = fig.colorbar(scatter, ax=axs[i, 1], pad=0.01)
        cbar.set_label('Phi error, °')
        
        axs[i, 1].set_xlim(-20, 20)
        axs[i, 1].set_ylim(-20, 20)
        axs[i, 1].set_aspect('equal')
        axs[i, 1].grid(True, linestyle='--', alpha=0.7)
        axs[i, 1].axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        axs[i, 1].axvline(x=0, color='gray', linestyle='-', alpha=0.5)

    plt.tight_layout()
    output_filename = f"testing_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()
