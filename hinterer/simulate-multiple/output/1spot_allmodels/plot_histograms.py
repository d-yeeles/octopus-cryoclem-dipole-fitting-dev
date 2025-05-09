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
   
    fig, axs = plt.subplots(3, 6, figsize=(40, 15))

    fig.suptitle(f'Hinterer Simulations, θ={inclination}°', fontsize=16, y=0.99)

    for i, dataset in enumerate(datasets):

        dataset_inc = dataset[abs(dataset['inc_tru'] - inclination) <= 5]

        if dataset_inc.empty:
          continue

        IQR_multiplier = 100

        data = dataset_inc["para_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data, ddof=1)
        se = std / np.sqrt(np.size(filtered_data))
        sns.histplot(filtered_data, ax=axs[i, 0], kde=True, fill=True, binwidth=1)
        axs[i, 0].axvline(mean, color=dred, linewidth=1)
        axs[i, 0].axvspan(mean - std, mean + std, color='red', alpha=0.1)
        axs[i, 0].set_xlim(-20, 20)
        axs[i, 0].set_xlabel('Δ$\parallel$, nm')
        axs[i, 0].set_ylabel('Counts')
        axs[i, 0].set_title(f"{model_names[i]}\nParallel localisation residuals\n" 
                            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

        data = dataset_inc["perp_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        sns.histplot(filtered_data, ax=axs[i, 1], kde=True, fill=True, binwidth=1)
        axs[i, 1].axvline(mean, color=dred, linewidth=1)
        axs[i, 1].axvspan(mean - std, mean + std, color='red', alpha=0.1)
        axs[i, 1].set_xlim(-20, 20)
        axs[i, 1].set_xlabel('Δ$\perp$, nm')
        axs[i, 1].set_ylabel('Counts')
        axs[i, 1].set_title(f"{model_names[i]}\nPerpendicular localisation residuals\n" 
                            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

        data = dataset_inc["inc_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        # Check if all data would fit in a single bin with binwidth=2.5
        data_range = max(data) - min(data)
        if data_range <= 2.5:  # If all data would fit in one bin with binwidth=2.5
            sns.histplot(data=data, ax=axs[i, 2], kde=True, fill=True, bins=1)  # Force 5 bins
        else:
            sns.histplot(data=data, ax=axs[i, 2], kde=True, fill=True, binwidth=2.5)
        axs[i, 2].axvline(mean, color=dred, linewidth=1)
        axs[i, 2].axvspan(mean - std, mean + std, color='red', alpha=0.1)
        axs[i, 2].set_xlim([-5, 95])
        axs[i, 2].set_xlabel('Δθ, °')
        axs[i, 2].set_ylabel('Counts')
        axs[i, 2].set_title(f"{model_names[i]}\nPolar residuals\n" 
                            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

        data = dataset_inc["az_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        # Check if all data would fit in a single bin with binwidth=2.5
        data_range = max(data) - min(data)
        if data_range <= 5:  # If all data would fit in one bin with binwidth=2.5
            sns.histplot(data=data, ax=axs[i, 3], kde=True, fill=True, bins=1)  # Force 5 bins
        else:
            sns.histplot(data=data, ax=axs[i, 3], kde=True, fill=True, binwidth=5)
        axs[i, 3].axvline(mean, color=dred, linewidth=1)
        axs[i, 3].axvspan(mean - std, mean + std, color='red', alpha=0.1)
        axs[i, 3].set_xlim([-10, 190])
        axs[i, 3].set_xlabel('Δ$\phi$, °')
        axs[i, 3].set_ylabel('Counts')
        axs[i, 3].set_title(f"{model_names[i]}\nAzimuth residuals\n" 
                            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

        data = dataset_inc["photon_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        # Check if all data would fit in a single bin with binwidth=2.5
        data_range = max(data) - min(data)
        if data_range <= 50:  # If all data would fit in one bin with binwidth=2.5
            sns.histplot(data=data, ax=axs[i, 4], kde=True, fill=True, bins=1)  # Force 5 bins
        else:
            sns.histplot(data=data, ax=axs[i, 4], kde=True, fill=True, binwidth=50)
        axs[i, 4].axvline(mean, color=dred, linewidth=1)
        axs[i, 4].axvspan(mean - std, mean + std, color='red', alpha=0.1)
        axs[i, 4].set_xlim([-220, 1220])
        axs[i, 4].set_xlabel('ΔN, °')
        axs[i, 4].set_ylabel('Counts')
        axs[i, 4].set_title(f"{model_names[i]}\nPhoton count residuals\n" 
                            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

        data = dataset_inc["obj_est"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        # Check if all data would fit in a single bin with binwidth=2.5
        data_range = max(data) - min(data)
        if data_range <= 100:  # If all data would fit in one bin with binwidth=2.5
            sns.histplot(data=data, ax=axs[i, 5], kde=True, fill=True, bins=1)  # Force 5 bins
        else:
            sns.histplot(data=data, ax=axs[i, 5], kde=True, fill=True, binwidth=100)
        axs[i, 5].axvline(mean, color=dred, linewidth=1)
        axs[i, 5].axvspan(mean - std, mean + std, color='red', alpha=0.1)
        axs[i, 5].set_xlim([0, 1200])
        axs[i, 5].set_xlabel('Log-likelihood, °')
        axs[i, 5].set_ylabel('Counts')
        axs[i, 5].set_title(f"{model_names[i]}\nLog-likelihood\n" 
                            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")


    plt.tight_layout()
    output_filename = f"histograms_hinterer_inc{round(inclination)}.pdf"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()

