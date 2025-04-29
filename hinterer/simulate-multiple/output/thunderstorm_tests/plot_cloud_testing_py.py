import numpy as np
import os
import sys
import importlib.util
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.image as mpimg
import subprocess
import seaborn as sns


def convert_lists_to_degrees(module, variable_names):
    for name in variable_names:
        if name in vars(module):
            vars(module)[name] = [x * 180 / np.pi for x in vars(module)[name]]

def compute_para_perp(x_diffs, y_diffs, azimuths):

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

# Initialise empty dataframes
data_gaussian = pd.DataFrame(columns=[
    "x_thu",
    "y_thu",
    "inc_thu",
    "az_thu",
    "x_est",
    "y_est",
    "inc_est",
    "az_est",
    "x_dif",
    "y_dif",
    "inc_dif",
    "az_dif",
    "para_dif",
    "perp_dif",
    "photon_thu",
    "photon_est",
    "photon_dif",
    "obj_est",
])

data_hinterer = data_gaussian.copy()
data_mortensen = data_gaussian.copy()

datasets = [data_gaussian, data_hinterer]
model_names = ['hinterer', 'mortensen']
file_paths = [
    './fitting_results_mortensen_test_hardcoded.py',
    './fitting_results_mortensen_test.py',
#    './fitting_results_mortensen_test_runtime.py',
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
                "x_thu": module.x_thu,
                "y_thu": module.y_thu,
                "inc_thu": module.inc_thu,
                "az_thu": module.az_thu,
                "x_est": module.x_est,
                "y_est": module.y_est,
                "inc_est": module.inc_est,
                "az_est": module.az_est,
                "x_dif": np.array(module.x_thu) - np.array(module.x_est),
                "y_dif": np.array(module.y_thu) - np.array(module.y_est),
                "inc_dif": np.array(module.inc_thu) - np.array(module.inc_est),
                "az_dif": np.array(module.az_thu) - np.array(module.az_est),
                "para_dif": compute_para_perp(np.array(module.x_thu) - np.array(module.x_est), np.array(module.y_thu) - np.array(module.y_est), module.az_thu)[0],
                "perp_dif": compute_para_perp(np.array(module.x_thu) - np.array(module.x_est), np.array(module.y_thu) - np.array(module.y_est), module.az_thu)[1],
                "photon_thu": module.photon_thu,
                "photon_est": module.photon_est,
                "photon_dif": module.photon_dif,
                "obj_est": module.obj_est,
            })

        if dataset.empty:
            datasets[i] = newdata
        else:
            datasets[i] = pd.concat([dataset, newdata], ignore_index=True)




for dataset in datasets:

    # Convert rad to deg
    angle_columns = ["inc_thu", "az_thu", "inc_est", "az_est", "inc_dif", "az_dif"]
    dataset[angle_columns] = dataset[angle_columns] * 180 / np.pi

    # Account for wrapping of inc around 180 degrees:
    dataset["inc_dif"] = np.mod(dataset["inc_dif"], 180)
    dataset["inc_dif"] = np.minimum(np.abs(dataset["inc_dif"]), 180 - np.abs(dataset["inc_dif"]))

    # Account for wrapping of az around 360 degrees:
    dataset["az_dif"] = np.mod(dataset["az_dif"], 360)
    dataset["az_dif"] = np.minimum(np.abs(dataset["az_dif"]), 360 - np.abs(dataset["az_dif"]))




output_dir = './'

import seaborn as sns

# Generate plots for each fixed inclination
for inclination in [68]:
   
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    for i, dataset in enumerate(datasets):

        dataset_inc = dataset[abs(dataset['inc_thu'] - inclination) <= 5]

        if dataset_inc.empty:
          continue

        IQR_multiplier = 500000000000
        
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
        
        # Filter all relevant data using the mask
        filtered_para = para_dif[mask]
        filtered_perp = perp_dif[mask]
        filtered_az_dif = dataset_inc["inc_dif"][mask]  # Using az_dif for coloring
        
        # Single scatter plot with para_dif on x-axis and perp_dif on y-axis
        scatter = axs[i, 0].scatter(filtered_para, filtered_perp, s=2, c=filtered_az_dif, cmap='rainbow')

#        # Add text labels
#        for j, (x, y, inc) in enumerate(zip(filtered_para, filtered_perp, filtered_az_dif)):
#            # Round inc_est to 1 decimal place
#            inc_label = f"{inc:.6f}"
#            axs[i, 0].annotate(inc_label, (x, y), fontsize=6, color='black')
 
        # Add axis labels and title
        axs[i, 0].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[i, 0].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[i, 0].set_title(
            f"{model_names[i]}\nθ={inclination}°"
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









        IQR_multiplier = 500000000000
        
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
        
        # Filter all relevant data using the mask
        filtered_para = para_dif[mask]
        filtered_perp = perp_dif[mask]
        filtered_az_dif = dataset_inc["az_dif"][mask]  # Using az_dif for coloring
        
        # Single scatter plot with para_dif on x-axis and perp_dif on y-axis
        scatter = axs[i, 1].scatter(filtered_para, filtered_perp, s=2, c=filtered_az_dif, cmap='rainbow')

#        # Add text labels
#        for j, (x, y, inc) in enumerate(zip(filtered_para, filtered_perp, filtered_az_dif)):
#            # Round inc_est to 1 decimal place
#            inc_label = f"{inc:.6f}"
#            axs[i, 1].annotate(inc_label, (x, y), fontsize=6, color='black')
 
        # Add axis labels and title
        axs[i, 1].set_xlabel('Parallel error (Δ$\parallel$), nm')
        axs[i, 1].set_ylabel('Perpendicular error (Δ$\perp$), nm')
        axs[i, 1].set_title(
            f"{model_names[i]}\nθ={inclination}°"
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
    output_filename = f"testing_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()

