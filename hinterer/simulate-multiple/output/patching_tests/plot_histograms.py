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
data_gaussian = pd.DataFrame(columns=[
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
    "photon_tru",
    "photon_est",
    "photon_err"
])

datasets = [data_gaussian]
model_names = ['hinterer 9spot']
file_paths = [
    './results_23spot_large.py',
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
                "photon_err": module.photon_err
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




output_dir = './results_23spot_large/'

import seaborn as sns

# Generate plots for each fixed inclination
inclinations = [45]#[0, 23, 45, 68, 90]
for inclination in inclinations:
   
    fig, axs = plt.subplots(1, 5, figsize=(35, 5))

    dataset_inc = datasets[0][abs(datasets[0]['inc_tru'] - inclination) <= 5]

    # Filter outliers
    iqr_multiplier = 5000000
    for column in ["para_err", "perp_err", "inc_err", "az_err", "photon_err"]:
        data = dataset_inc[column]
        q1 = np.percentile(data, 25)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        lower_bound = q1 - iqr_multiplier * iqr
        upper_bound = q3 + iqr_multiplier * iqr
        dataset_inc[column] = data[(data >= lower_bound) & (data <= upper_bound)]

    data = dataset_inc["para_err"]
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    se = std / np.sqrt(np.size(data))
    sns.histplot(data, ax=axs[0], kde=True, fill=True, binwidth=2)
    axs[0].axvline(mean, color=dred, linewidth=1)
    axs[0].axvspan(mean - std, mean + std, color='red', alpha=0.1)
    axs[0].set_xlim(-45, 45)
    axs[0].set_xlabel('Δ$\parallel$, nm')
    axs[0].set_ylabel('Counts')
    axs[0].set_title(f"{model_names[0]}\nParallel localisation residuals\n"
        f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

    data = dataset_inc["perp_err"]
    mean = np.mean(data)
    std = np.std(data)
    se = std / np.sqrt(np.size(data))
    sns.histplot(data, ax=axs[1], kde=True, fill=True, binwidth=2)
    axs[1].axvline(mean, color=dred, linewidth=1)
    axs[1].axvspan(mean - std, mean + std, color='red', alpha=0.1)
    axs[1].set_xlim(-40, 40)
    axs[1].set_xlabel('Δ$\perp$, nm')
    axs[1].set_ylabel('Counts')
    axs[1].set_title(f"{model_names[i]}\nPerpendicular localisation residuals\n" 
                        f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

    data = dataset_inc["inc_err"]
    mean = np.mean(data)
    std = np.std(data)
    se = std / np.sqrt(np.size(data))
    sns.histplot(data=data, ax=axs[2], kde=True, fill=True, binwidth=2.5)
    axs[2].axvline(mean, color=dred, linewidth=1)
    axs[2].axvspan(mean - std, mean + std, color='red', alpha=0.1)
    axs[2].set_xlim([-5, 95])
    axs[2].set_xlabel('Δθ, °')
    axs[2].set_ylabel('Counts')
    axs[2].set_title(f"{model_names[i]}\nPolar residuals\n" 
                        f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

    data = dataset_inc["az_err"]
    mean = np.mean(data)
    std = np.std(data)
    se = std / np.sqrt(np.size(data))
    sns.histplot(data=data, ax=axs[3], kde=True, fill=True, binwidth=5)
    axs[3].axvline(mean, color=dred, linewidth=1)
    axs[3].axvspan(mean - std, mean + std, color='red', alpha=0.1)
    axs[3].set_xlim([-10, 190])
    axs[3].set_xlabel('Δ$\phi$, °')
    axs[3].set_ylabel('Counts')
    axs[3].set_title(f"{model_names[i]}\nAzimuth residuals\n" 
                        f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

    data = dataset_inc["photon_err"]
    mean = np.mean(data)
    std = np.std(data)
    se = std / np.sqrt(np.size(data))
    sns.histplot(data=data, ax=axs[4], kde=True, fill=True, binwidth=50)
    axs[4].axvline(mean, color=dred, linewidth=1)
    axs[4].axvspan(mean - std, mean + std, color='red', alpha=0.1)
    axs[4].set_xlim([-220, 1220])
    axs[4].set_xlabel('ΔN, °')
    axs[4].set_ylabel('Counts')
    axs[4].set_title(f"{model_names[i]}\nPhoton count residuals\n" 
                        f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}")

    plt.tight_layout()
    output_filename = f"histograms_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()

