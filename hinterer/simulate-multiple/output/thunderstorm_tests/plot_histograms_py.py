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

datasets = [data_gaussian, data_hinterer, data_mortensen]
model_names = ['gaussian on mortensen', 'hinterer on mortensen', 'mortensen on mortensen']
file_paths = [
#    './fitting_results_gaussian_on_hinterer_all.py',
    './fitting_results_hinterer_on_mortensen_all.py',
    './fitting_results_hinterer_on_mortensen_all.py',
    './fitting_results_mortensen_on_mortensen_all.py',
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
for inclination in [0]:
   
    fig, axs = plt.subplots(3, 6, figsize=(40, 15))

    for i, dataset in enumerate(datasets):

        dataset_inc = dataset[abs(dataset['inc_thu'] - inclination) <= 5]

        if dataset_inc.empty:
          continue

        IQR_multiplier = 500000000

        data = dataset_inc["para_dif"]
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

        data = dataset_inc["perp_dif"]
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

        data = dataset_inc["inc_dif"]
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
        sns.histplot(data=data, ax=axs[i, 2], kde=True, fill=True, binwidth=2.5)
        axs[i, 2].axvline(mean, color=dred, linewidth=1)
        axs[i, 2].axvspan(mean - std, mean + std, color='red', alpha=0.1)
        axs[i, 2].set_xlim([-5, 95])
        axs[i, 2].set_xlabel('Δθ, °')
        axs[i, 2].set_ylabel('Counts')
        axs[i, 2].set_title(f"{model_names[i]}\nPolar residuals\n" 
                            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

        data = dataset_inc["az_dif"]
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
        sns.histplot(data=data, ax=axs[i, 3], kde=True, fill=True, binwidth=5)
        axs[i, 3].axvline(mean, color=dred, linewidth=1)
        axs[i, 3].axvspan(mean - std, mean + std, color='red', alpha=0.1)
        axs[i, 3].set_xlim([-10, 190])
        axs[i, 3].set_xlabel('Δ$\phi$, °')
        axs[i, 3].set_ylabel('Counts')
        axs[i, 3].set_title(f"{model_names[i]}\nAzimuth residuals\n" 
                            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

        data = dataset_inc["photon_dif"]
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
        sns.histplot(data=data, ax=axs[i, 5], kde=True, fill=True, binwidth=100)
        axs[i, 5].axvline(mean, color=dred, linewidth=1)
        axs[i, 5].axvspan(mean - std, mean + std, color='red', alpha=0.1)
        axs[i, 5].set_xlim([0, 15000])
        axs[i, 5].set_xlabel('Log-likelihood, °')
        axs[i, 5].set_ylabel('Counts')
        axs[i, 5].set_title(f"{model_names[i]}\nLog-likelihood\n" 
                            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")


    plt.tight_layout()
    output_filename = f"histograms_hinterer_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()

