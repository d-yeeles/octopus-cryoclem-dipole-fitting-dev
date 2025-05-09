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
    "photon_err",
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
                "photon_err": module.photon_err,
                "obj_est": module.obj_est,
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

import seaborn as sns

# Generate plots for each fixed inclination
for inclination in [0]:
   
    fig, axs = plt.subplots(3, 6, figsize=(40, 15))

    for i, dataset in enumerate(datasets):

        dataset_inc = dataset[abs(dataset['inc_tru'] - inclination) <= 5]

        if dataset_inc.empty:
          continue

        IQR_multiplier = 500000000

        data = dataset_inc["para_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_az = dataset_inc["az_tru"][(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data, ddof=1)
        se = std / np.sqrt(np.size(filtered_data))
        axs[i, 0].scatter(filtered_az, filtered_data, s=2)
        axs[i, 0].axhline(mean, color=dred, linewidth=1)
        axs[i, 0].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 0].set_ylim(-20, 20)
        axs[i, 0].set_xlabel('$\phi$, °')
        axs[i, 0].set_ylabel('Δ$\parallel$, nm')
        axs[i, 0].set_title(
            f"{model_names[i]}\nParallel localisation residuals\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )
        #axs[i, 0].legend(loc='upper right')

        data = dataset_inc["perp_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_az = dataset_inc["az_tru"][(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        axs[i, 1].scatter(filtered_az, filtered_data, s=2)
        axs[i, 1].axhline(mean, color=dred, linewidth=1)#, label=f'μ = {mean:.2f}\nσ = {std:.2f}')
        axs[i, 1].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 1].set_ylim(-20, 20)
        axs[i, 1].set_xlabel('$\phi$, °')
        axs[i, 1].set_ylabel('Δ$\perp$, nm')
        axs[i, 1].set_title(
            f"{model_names[i]}\nPerpendicular localisation residuals\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )
        #axs[i, 1].legend(loc='upper right')

        data = dataset_inc["inc_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_az = dataset_inc["az_tru"][(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        axs[i, 2].scatter(filtered_az, filtered_data, s=2)
        axs[i, 2].axhline(mean, color=dred, linewidth=1)#, label=f'μ = {mean:.2f}\nσ = {std:.2f}')
        axs[i, 2].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 2].set_ylim(0, 90)
        axs[i, 2].set_xlabel('$\phi$, °')
        axs[i, 2].set_ylabel('Δθ, °')
        #axs[i, 2].set_title(model_names[i]+"\nPolar residuals")
        axs[i, 2].set_title(
            f"{model_names[i]}\nPolar residuals\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )
        #axs[i, 2].legend(loc='upper right')

        data = dataset_inc["az_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_az = dataset_inc["az_tru"][(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        axs[i, 3].scatter(filtered_az, filtered_data, s=2)
        axs[i, 3].axhline(mean, color=dred, linewidth=1)#, label=f'μ = {mean:.2f}\nσ = {std:.2f}')
        axs[i, 3].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 3].set_ylim(0, 180)
        axs[i, 3].set_xlabel('$\phi$, °')
        axs[i, 3].set_ylabel('Δ$\phi$, °')
        axs[i, 3].set_title(
            f"{model_names[i]}\nAzimuth residuals\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )
        #axs[i, 3].legend(loc='upper right')

        data = dataset_inc["photon_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_az = dataset_inc["az_tru"][(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        axs[i, 4].scatter(filtered_az, filtered_data, s=2)
        axs[i, 4].axhline(mean, color=dred, linewidth=1)#, label=f'μ = {mean:.2f}\nσ = {std:.2f}')
        axs[i, 4].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 4].set_ylim(-200, 1200)
        axs[i, 4].set_xlabel('$\phi$, °')
        axs[i, 4].set_ylabel('ΔN')
        axs[i, 4].set_title(
            f"{model_names[i]}\nPhoton count residuals\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )

        data = dataset_inc["obj_est"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_az = dataset_inc["az_tru"][(data >= lower_bound) & (data <= upper_bound)]
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        axs[i, 5].scatter(filtered_az, filtered_data, s=2)
        axs[i, 5].axhline(mean, color=dred, linewidth=1)#, label=f'μ = {mean:.2f}\nσ = {std:.2f}')
        axs[i, 5].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 5].set_ylim(0, 15000)
        axs[i, 5].set_xlabel('$\phi$, °')
        axs[i, 5].set_ylabel('log-likelihood')
        axs[i, 5].set_title(
            f"{model_names[i]}\nLog-likelihood\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )

    plt.tight_layout()
    output_filename = f"scatter_hinterer_inc{round(inclination)}.pdf"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()

