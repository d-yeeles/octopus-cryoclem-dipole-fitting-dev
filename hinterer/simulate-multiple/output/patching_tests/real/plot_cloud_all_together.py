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
data_1 = pd.DataFrame(columns=[
    "x_tru",
    "y_tru",
#    "inc_tru",
#    "az_tru",
    "x_est",
    "y_est",
    "inc_est",
    "az_est",
    "x_err",
    "y_err",
#    "inc_err",
#    "az_err",
#    "para_err",
#    "perp_err",
#    "obj_tru",
#    "obj_est",
#    "obj_err"
#    "photon_tru",
    "photon_est",
#    "photon_err"
])

datasets = [data_1]
model_names = ['hinterer vs thunderstorm']
file_paths = [
    './results_real.py',
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
#                "inc_tru": module.inc_tru,
#                "az_tru": module.az_tru,
                "x_est": module.x_est,
                "y_est": module.y_est,
                "inc_est": module.inc_est,
                "az_est": module.az_est,
                "x_err": np.array(module.x_tru) - np.array(module.x_est),
                "y_err": np.array(module.y_tru) - np.array(module.y_est),
#                "inc_err": np.array(module.inc_tru) - np.array(module.inc_est),
#                "az_err": np.array(module.az_tru) - np.array(module.az_est),
#                "para_err": compute_para_perp(np.array(module.x_tru) - np.array(module.x_est), np.array(module.y_tru) - np.array(module.y_est), module.az_tru)[0],
#                "perp_err": compute_para_perp(np.array(module.x_tru) - np.array(module.x_est), np.array(module.y_tru) - np.array(module.y_est), module.az_tru)[1],
#                "obj_tru": module.obj_tru,
#                "obj_est": module.obj_est,
#                "obj_err": module.obj_err
#                "photon_tru": module.photon_tru,
                "photon_est": module.photon_est,
#                "photon_err": np.array(module.photon_tru) - np.array(module.photon_est)

            })

        if dataset.empty:
            datasets[i] = newdata
        else:
            datasets[i] = pd.concat([dataset, newdata], ignore_index=True)

for dataset in datasets:

    # Convert rad to deg
    angle_columns = ["inc_est", "az_est"]
    dataset[angle_columns] = dataset[angle_columns] * 180 / np.pi

output_dir = './'

# Plots

fig, axs = plt.subplots(1, 1, figsize=(12, 12))


dataset_inc = datasets[0]

scatter = axs.scatter(dataset_inc["x_err"], dataset_inc["y_err"], s=2)
 
#axs.scatter(dataset_inc["x_tru"], dataset_inc["y_tru"], s=20, color='red', label='True positions')
#axs.scatter(dataset_inc["x_est"], dataset_inc["y_est"], s=20, color='green', label='Estimated positions')

## Draw black lines connecting true and estimated positions
#for i in range(len(dataset_inc)):
#    axs.plot([dataset_inc["x_tru"].iloc[i], dataset_inc["x_est"].iloc[i]], 
#             [dataset_inc["y_tru"].iloc[i], dataset_inc["y_est"].iloc[i]], 
#             'k-', linewidth=0.5)
   
# Add axis labels and title
axs.set_xlabel('Δx, nm')
axs.set_ylabel('Δy, nm')
axs.set_title(
    f"{model_names[0]}°\n"
)
    
axs.set_xlim(-100, 100)#(-max_val*1.1, max_val*1.1)
axs.set_ylim(-100, 100)#(-max_val*1.1, max_val*1.1)
    
# Force aspect ratio to be 1:1
axs.set_aspect('equal')
axs.grid(True, linestyle='--', alpha=0.7)
axs.set_xticks(np.arange(-100, 100, 20))
axs.set_yticks(np.arange(-100, 100, 20))
axs.axhline(y=0, color='gray', linestyle='-', alpha=0.5)
axs.axvline(x=0, color='gray', linestyle='-', alpha=0.5)
    
plt.tight_layout()
output_filename = f"visualisation.png"
plt.savefig(f"{output_dir}{output_filename}", dpi=300)
plt.close()
