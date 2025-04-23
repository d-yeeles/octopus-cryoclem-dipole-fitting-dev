import numpy as np
import os
import sys
import importlib.util
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.image as mpimg
import subprocess


def convert_lists_to_degrees(module, variable_names):
    for name in variable_names:
        if name in vars(module):
            vars(module)[name] = [x * 180 / np.pi for x in vars(module)[name]]

def compute_para_perp(x_errors, y_errors, azimuths):
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

# Initialise empty dataframes
data_mortensen = pd.DataFrame(columns=[
    "x_tru", "y_tru", "inc_tru", "az_tru", "x_est", "y_est", "inc_est", "az_est",
    "x_err", "y_err", "inc_err", "az_err", "para_err", "perp_err",
    "photon_tru", "photon_est", "photon_err"
])

data_2 = data_mortensen.copy()
data_3 = data_mortensen.copy()
data_4 = data_mortensen.copy()
data_5 = data_mortensen.copy()
data_6 = data_mortensen.copy()

datasets = [data_mortensen, data_2, data_3, data_4, data_5, data_6]
model_names = [
    'Gaussian fit on Hinterer sims',
    'Hinterer fit on Hinterer sims',
    'Mortensen fit on Hinterer sims',
    'Gaussian fit on Mortensen sims',
    'Hinterer fit on Mortensen sims',
    'Mortensen fit on Mortensen sims',
]
file_paths = [
    'results_gaussian_on_hinterer.py',
    'results_hinterer_on_hinterer.py',
    'results_mortensen_on_hinterer.py',
    'results_gaussian_on_mortensen.py',
    'results_hinterer_on_mortensen.py',
    'results_mortensen_on_mortensen.py',
]


# Importing all the data from each results file
subdirectory = "./"

for i, dataset in enumerate(datasets):
    file_path = os.path.join(subdirectory, file_paths[i])

    if os.path.exists(file_path):
        spec = importlib.util.spec_from_file_location(model_names[i], file_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        newdata = pd.DataFrame({
            "x_tru": module.x_tru,
            "y_tru": module.y_tru,
            "inc_tru": module.inc_tru,
            "az_tru": module.az_tru,
            "x_est": module.x_est,
            "y_est": module.y_est,
            "inc_est": module.inc_est,
            "az_est": module.az_est,
            "x_err": np.array(module.x_tru) - np.array(module.x_est),
            "y_err": np.array(module.y_tru) - np.array(module.y_est),
            "inc_err": np.array(module.inc_tru) - np.array(module.inc_est),
            "az_err": np.array(module.az_tru) - np.array(module.az_est),
            "para_err": compute_para_perp(np.array(module.x_tru) - np.array(module.x_est), np.array(module.y_tru) - np.array(module.y_est), module.az_tru)[0],
            "perp_err": compute_para_perp(np.array(module.x_tru) - np.array(module.x_est), np.array(module.y_tru) - np.array(module.y_est), module.az_tru)[1],
            "photon_tru": module.photon_tru,
            "photon_est": module.photon_est,
            "photon_err": np.array(module.photon_tru) - np.array(module.photon_est)
        })

        if dataset.empty:
            datasets[i] = newdata
        else:
            datasets[i] = pd.concat([dataset, newdata], ignore_index=True)


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

# Reorganize datasets by simulation type
hinterer_datasets = []  # 'on Hinterer' datasets for top row
mortensen_datasets = []  # 'on Mortensen' datasets for bottom row

for i, name in enumerate(model_names):
    if 'on Hinterer' in name:
        hinterer_datasets.append((datasets[i], name))
    elif 'on Mortensen' in name:
        mortensen_datasets.append((datasets[i], name))

# Generate plots for each fixed inclination
for i, inclination in enumerate([0, 23, 45, 68, 90]):

    # TRANSPOSED: Change from (2, 3) to (3, 2)
    fig, axs = plt.subplots(3, 2, figsize=(12, 18))
    plt.suptitle(f'θ = {inclination}', fontsize=16, y=0.98)

    # Plot the 'on Hinterer' datasets in the left column (was top row)
    for j, (dataset, name) in enumerate(hinterer_datasets):
        dataset_inc = dataset[np.abs(dataset['inc_tru'] - inclination) <= 5]

        # TRANSPOSED: Change indexing from [0, j] to [j, 0]
        scatter = axs[j, 0].scatter(dataset_inc["para_err"], dataset_inc["perp_err"], s=2, c=dataset_inc["az_err"], cmap='rainbow')
        
        # Add axis labels and title
        axs[j, 0].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[j, 0].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[j, 0].set_title(name)
        
        # Add colorbar
        cbar = fig.colorbar(scatter, ax=axs[j, 0], pad=0.01)
        cbar.set_label('Phi error')
        
        axs[j, 0].set_xlim(-50, 50)
        axs[j, 0].set_ylim(-50, 50)
        
        # Force aspect ratio to be 1:1
        axs[j, 0].set_aspect('equal')
        axs[j, 0].grid(True, linestyle='--', alpha=0.7)
        axs[j, 0].axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        axs[j, 0].axvline(x=0, color='gray', linestyle='-', alpha=0.5)

    # Plot the 'on Mortensen' datasets in the right column (was bottom row)
    for j, (dataset, name) in enumerate(mortensen_datasets):
        dataset_inc = dataset[np.abs(dataset['inc_tru'] - inclination) <= 5]

        # TRANSPOSED: Change indexing from [1, j] to [j, 1]
        scatter = axs[j, 1].scatter(dataset_inc["para_err"], dataset_inc["perp_err"], s=2, c=dataset_inc["az_err"], cmap='rainbow')
        
        # Add axis labels and title
        axs[j, 1].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[j, 1].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[j, 1].set_title(name)
        
        # Add colorbar
        cbar = fig.colorbar(scatter, ax=axs[j, 1], pad=0.01)
        cbar.set_label('Phi error')
        
        axs[j, 1].set_xlim(-50, 50)
        axs[j, 1].set_ylim(-50, 50)
        
        # Force aspect ratio to be 1:1
        axs[j, 1].set_aspect('equal')
        axs[j, 1].grid(True, linestyle='--', alpha=0.7)
        axs[j, 1].axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        axs[j, 1].axvline(x=0, color='gray', linestyle='-', alpha=0.5)

    # TRANSPOSED: Handle case if there's an uneven number of plots in each column
    max_plots = max(len(hinterer_datasets), len(mortensen_datasets))
    for col in range(2):
        for row in range(len(hinterer_datasets if col == 0 else mortensen_datasets), 3):
            axs[row, col].set_visible(False)

    plt.tight_layout()
    output_filename = f"cloud_results_inc{inclination}.png"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()
