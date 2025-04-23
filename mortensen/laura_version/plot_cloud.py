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
    "x_tru",
    "y_tru",
    "inc_tru",
    "az_tru",
    "x_est",
    "y_est",
    "inc_est",
    "az_est",
    "x_err",
    "y_err",
    "inc_err",
    "az_err",
    "para_err",
    "perp_err",
#    "obj_tru",
#    "obj_est",
#    "obj_err"
    "photon_tru",
    "photon_est",
    "photon_err"
])

data_mortensen_2 = data_mortensen.copy()

datasets = [data_mortensen, data_mortensen_2]
model_names = ['mortensen', 'mortensen 2']
module_names = ['mortensen', 'mortensen 2']
file_paths = [
    'fitting_results_mortensen_4.py',
    'fitting_results_mortensen_5.py',
]



# Importing all the data from each results file
subdirectory = "./"

for i, dataset in enumerate(datasets):

    file_path = os.path.join(subdirectory, file_paths[i])

    if os.path.exists(file_path):

        spec = importlib.util.spec_from_file_location(module_names[i], file_path)
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
#                "obj_tru": module.obj_tru,
#                "obj_est": module.obj_est,
#                "obj_err": module.obj_err
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

# Azimuth plots

fig, axs = plt.subplots(2, 2, figsize=(12, 12))

# Generate plots for each fixed inclination
for i, inclination in enumerate([45, 68]):

    for j, dataset in enumerate(datasets):

        dataset_inc = dataset[np.abs(dataset['inc_tru'] - inclination) <= 5]

        IQR_multiplier = 5000
        
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
        scatter = axs[i, j].scatter(filtered_para, filtered_perp, s=2, c=filtered_az_err, cmap='rainbow')#, vmin=0, vmax=90)
    
#        # Add text labels for az_est beside each point
#        for x, y, az in zip(filtered_para, filtered_perp, filtered_az_err):
#            axs[i, j].text(x, y, f'{az:.1f}', fontsize=6, color='black', ha='left', va='bottom', alpha=0.7)

    
        # Add axis labels and title
        axs[i, j].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[i, j].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[i, j].set_title(
            f"{model_names[j]} \n θ = {inclination}°\n"
        )
        
        # Add colorbar
        cbar = fig.colorbar(scatter, ax=axs[i, j], pad=0.01)
        cbar.set_label('Phi err')
        
        # Add a reference circle showing the IQR*5 boundary
        max_val = max(abs(filtered_para.max()), abs(filtered_para.min()), 
                      abs(filtered_perp.max()), abs(filtered_perp.min()))
        axs[i, j].set_xlim(-max_val*1.1, max_val*1.1)
        axs[i, j].set_ylim(-max_val*1.1, max_val*1.1)
        
        # Force aspect ratio to be 1:1
        axs[i, j].set_aspect('equal')
        axs[i, j].grid(True, linestyle='--', alpha=0.7)
        axs[i, j].axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        axs[i, j].axvline(x=0, color='gray', linestyle='-', alpha=0.5)
   



    plt.tight_layout()
    output_filename = f"cloud_results_allinc.png"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()
