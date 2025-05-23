import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns


def compute_para_perp(x_errfss, y_errfss, azimuths):
    """Compute parallel and perpendicular errors"""
    para_errfss = x_errfss*np.cos(azimuths) - y_errfss*np.sin(azimuths)
    perp_errfss = x_errfss*np.sin(azimuths) + y_errfss*np.cos(azimuths)
    return para_errfss, perp_errfss


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
    'results_posterior_test.csv',
]

model_names = ['Hinterer']
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
            para_errfss, perp_errfss = compute_para_perp(
                df['x_err'],# - df['x_est'], 
                df['y_err'],# - df['y_est'], 
                az_rad
            )
            
            df['para_err'] = para_errfss
            df['perp_err'] = perp_errfss
        
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

fig, axs = plt.subplots(1, 6, figsize=(40, 5))

fig.suptitle(f'Mortensen Simulations', fontsize=16, y=0.99)

for i, dataset in enumerate(datasets):

    dataset_inc = dataset

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
    sns.histplot(filtered_data, ax=axs[0], kde=True, fill=True, binwidth=1)
    axs[0].axvline(mean, color=dred, linewidth=1)
    axs[0].axvspan(mean - std, mean + std, color='red', alpha=0.1)
    axs[0].set_xlim(-20, 20)
    axs[0].set_xlabel('Δ$\parallel$, nm')
    axs[0].set_ylabel('Counts')
    axs[0].set_title(f"{model_names[i]}\nParallel localisation residuals\n" 
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
    sns.histplot(filtered_data, ax=axs[1], kde=True, fill=True, binwidth=1)
    axs[1].axvline(mean, color=dred, linewidth=1)
    axs[1].axvspan(mean - std, mean + std, color='red', alpha=0.1)
    axs[1].set_xlim(-20, 20)
    axs[1].set_xlabel('Δ$\perp$, nm')
    axs[1].set_ylabel('Counts')
    axs[1].set_title(f"{model_names[i]}\nPerpendicular localisation residuals\n" 
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
        sns.histplot(data=data, ax=axs[2], kde=True, fill=True, bins=1)  # Force 5 bins
    else:
        sns.histplot(data=data, ax=axs[2], kde=True, fill=True, binwidth=2.5)
    axs[2].axvline(mean, color=dred, linewidth=1)
    axs[2].axvspan(mean - std, mean + std, color='red', alpha=0.1)
    axs[2].set_xlim([-5, 95])
    axs[2].set_xlabel('Δθ, °')
    axs[2].set_ylabel('Counts')
    axs[2].set_title(f"{model_names[i]}\nPolar residuals\n" 
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
        sns.histplot(data=data, ax=axs[3], kde=True, fill=True, bins=1)  # Force 5 bins
    else:
        sns.histplot(data=data, ax=axs[3], kde=True, fill=True, binwidth=5)
    axs[3].axvline(mean, color=dred, linewidth=1)
    axs[3].axvspan(mean - std, mean + std, color='red', alpha=0.1)
    axs[3].set_xlim([-10, 190])
    axs[3].set_xlabel('Δ$\phi$, °')
    axs[3].set_ylabel('Counts')
    axs[3].set_title(f"{model_names[i]}\nAzimuth residuals\n" 
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
    sns.histplot(data=data, ax=axs[4], kde=True, fill=True)#, binwidth=50)
    axs[4].axvline(mean, color=dred, linewidth=1)
    axs[4].axvspan(mean - std, mean + std, color='red', alpha=0.1)
#        axs[4].set_xlim([-220, 1220])
    axs[4].set_xlabel('ΔN, °')
    axs[4].set_ylabel('Counts')
    axs[4].set_title(f"{model_names[i]}\nPhoton count residuals\n" 
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
    sns.histplot(data=data, ax=axs[5], kde=True, fill=True)#, binwidth=100)
    axs[5].axvline(mean, color=dred, linewidth=1)
    axs[5].axvspan(mean - std, mean + std, color='red', alpha=0.1)
#        axs[5].set_xlim([0, 1200])
    axs[5].set_xlabel('Log-likelihood, °')
    axs[5].set_ylabel('Counts')
    axs[5].set_title(f"{model_names[i]}\nLog-likelihood\n" 
                        f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")


plt.tight_layout()
output_filename = f"testing_hist.png"
plt.savefig(f"{output_dir}{output_filename}", dpi=300)
plt.close()

