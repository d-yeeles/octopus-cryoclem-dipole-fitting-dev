import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns


def compute_para_perp(x_errors, y_errors, azimuths):
    """Compute parallel and perpendicular errors"""
    para_errors = x_errors*np.cos(azimuths) - y_errors*np.sin(azimuths)
    perp_errors = x_errors*np.sin(azimuths) + y_errors*np.cos(azimuths)
    return para_errors, perp_errors


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
    './fitting_results_gaussian_on_hinterer.csv',
    './fitting_results_hinterer_on_hinterer.csv',
    './fitting_results_mortensen_on_hinterer.csv',
]

model_names = ['Gaussian', 'Hinterer', 'Mortensen']
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
        
        datasets.append(df)
    else:
        print(f"Warning: File not found: {full_path}")
        # Add an empty DataFrame to maintain index alignment
        datasets.append(pd.DataFrame())

for dataset in datasets:

    # Convert rad to deg
    angle_columns = ["inc_tru", "az_tru", "inc_est", "az_est", "inc_err", "az_err"]
    dataset[angle_columns] = dataset[angle_columns] * 180 / np.pi

    # Account for wrapping of inc around 180 degrees:
    dataset["inc_err"] = np.mod(dataset["inc_err"], 180)
    dataset["inc_err"] = np.minimum(np.abs(dataset["inc_err"]), 180 - np.abs(dataset["inc_err"]))

    # Account for wrapping of az around 360 degrees:
    dataset["az_err"] = np.mod(dataset["az_err"], 360)
    dataset["az_err"] = np.minimum(np.abs(dataset["az_err"]), 360 - np.abs(dataset["az_err"]))




output_dir = './'

import seaborn as sns

# Generate plots for each fixed inclination
for inclination in [0, 23, 45, 68, 90]:
   
    fig, axs = plt.subplots(3, 2, figsize=(10, 15))

    fig.suptitle(f'Hinterer Simulations, θ={inclination}°', fontsize=16, y=0.99)

    for i, dataset in enumerate(datasets):

        dataset_inc = dataset[abs(dataset['inc_tru'] - inclination) <= 5]

        if dataset_inc.empty:
          continue

        IQR_multiplier = 100
        
        # Get error data
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
        filtered_az_err = dataset_inc["inc_err"][mask]  # Using az_err for coloring
        
        # Single scatter plot with para_err on x-axis and perp_err on y-axis
        scatter = axs[i, 0].scatter(filtered_para, filtered_perp, s=2, c=filtered_az_err, cmap='rainbow')

#        # Add text labels
#        for j, (x, y, inc) in enumerate(zip(filtered_para, filtered_perp, filtered_az_err)):
#            # Round inc_est to 1 decimal place
#            inc_label = f"{inc:.6f}"
#            axs[i, 0].annotate(inc_label, (x, y), fontsize=6, color='black')
 
        # Add axis labels and title
        axs[i, 0].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[i, 0].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[i, 0].set_title(
            f"{model_names[i]}"
        )
        
        # Add colorbar
        cbar = fig.colorbar(scatter, ax=axs[i, 0], pad=0.01)
        cbar.set_label('Theta error')
        
        axs[i, 0].set_xlim(-20, 20)
        axs[i, 0].set_ylim(-20, 20)
        
        # Force aspect ratio to be 1:1
        axs[i, 0].set_aspect('equal')
        axs[i, 0].grid(True, linestyle='--', alpha=0.7)
        axs[i, 0].axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        axs[i, 0].axvline(x=0, color='gray', linestyle='-', alpha=0.5)









        IQR_multiplier = 100
        
        # Get error data
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
        filtered_az_err = dataset_inc["az_err"][mask]  # Using az_err for coloring
        
        # Single scatter plot with para_err on x-axis and perp_err on y-axis
        scatter = axs[i, 1].scatter(filtered_para, filtered_perp, s=2, c=filtered_az_err, cmap='rainbow')

#        # Add text labels
#        for j, (x, y, inc) in enumerate(zip(filtered_para, filtered_perp, filtered_az_err)):
#            # Round inc_est to 1 decimal place
#            inc_label = f"{inc:.6f}"
#            axs[i, 1].annotate(inc_label, (x, y), fontsize=6, color='black')
 
        # Add axis labels and title
        axs[i, 1].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[i, 1].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[i, 1].set_title(
            f"{model_names[i]}"
        )
        
        # Add colorbar
        cbar = fig.colorbar(scatter, ax=axs[i, 1], pad=0.01)
        cbar.set_label('Phi error')
        
        axs[i, 1].set_xlim(-20, 20)
        axs[i, 1].set_ylim(-20, 20)
        
        # Force aspect ratio to be 1:1
        axs[i, 1].set_aspect('equal')
        axs[i, 1].grid(True, linestyle='--', alpha=0.7)
        axs[i, 1].axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        axs[i, 1].axvline(x=0, color='gray', linestyle='-', alpha=0.5)






    plt.tight_layout()
    output_filename = f"cloud_hinterer_inc{round(inclination)}.pdf"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()

