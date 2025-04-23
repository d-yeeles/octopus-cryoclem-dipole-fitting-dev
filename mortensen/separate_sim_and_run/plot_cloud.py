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
model_names = ['mortensen together', 'mortensen split']
module_names = ['mortensen', 'mortensen 2']
file_paths = [
#    'fitting_results_mortensen_allinone.py',
#    'fitting_results_mortensen_split_joined.py',
#    'fitting_results_mortensen_split.py',
#    'fitting_results_mortensen_array_2.py',
    'fitting_results_mortensen_on_itself.py',
    'fitting_results_mortensen_on_itself.py',
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

# Generate plots for each fixed inclination
for i, inclination in enumerate([0, 23, 45, 68, 90]):

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    for j, dataset in enumerate(datasets):

        dataset_inc = dataset[np.abs(dataset['inc_tru'] - inclination) <= 5]

        # Single scatter plot with para_err on x-axis and perp_err on y-axis
        scatter = axs[j].scatter(dataset_inc["para_err"], dataset_inc["perp_err"], s=2, c=dataset_inc["inc_est"], cmap='rainbow') # Using phi err colouring
    
#        # Add text labels for az_est beside each poinit
#        for x, y, az in zip(dataset_inc["para_err"], dataset_inc["perp_err"], dataset_inc["az_err"]):
#            axs[j].text(x, y, f'{az:.1f}', fontsize=6, color='black', ha='left', va='bottom', alpha=0.7)

        # Add axis labels and title
        axs[j].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[j].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[j].set_title(
            f"{model_names[j]} \n θ = {inclination}°\n"
        )
        
        # Add colorbar
        cbar = fig.colorbar(scatter, ax=axs[j], pad=0.01)
        cbar.set_label('Theta estimate')
        
#        max_val = max(abs(dataset_inc["para_err"].max()), abs(dataset_inc["para_err"].min()), 
#                      abs(dataset_inc["perp_err"].max()), abs(dataset_inc["perp_err"].min()))
        axs[j].set_xlim(-30,30)#-max_val*1.1, max_val*1.1)
        axs[j].set_ylim(-30,30)#-max_val*1.1, max_val*1.1)
        
        # Force aspect ratio to be 1:1
        axs[j].set_aspect('equal')
        axs[j].grid(True, linestyle='--', alpha=0.7)
        axs[j].axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        axs[j].axvline(x=0, color='gray', linestyle='-', alpha=0.5)
   



    plt.tight_layout()
    output_filename = f"cloud_results_inc{inclination}.png"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()
