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
data_gaussian_fm = pd.DataFrame(columns=[
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
    "obj_tru",
    "obj_est",
    "obj_err"
])

data_gaussian_ps = data_gaussian_fm.copy()
data_hinterer_fm = data_gaussian_fm.copy()
data_hinterer_ps = data_gaussian_fm.copy()

datasets = [data_gaussian_fm, data_gaussian_ps, data_hinterer_fm, data_hinterer_ps]
model_names = ['nothing', 'nothing', 'mask', 'patch']
module_names = ['gaussian_fm', 'gaussian_ps', 'hinterer_fm', 'hinterer_ps']
file_paths = [
    'fitting_results.py',
    '../1spot_vary_az_reduced/fitting_results_blah.py',
    'fitting_results_mask.py',
    'fitting_results_patch.py',
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
                "x_err": module.x_err,
                "y_err": module.y_err,
                "inc_err": module.inc_err,
                "az_err": module.az_err,
                "para_err": compute_para_perp(module.x_err, module.y_err, module.az_tru)[0],
                "perp_err": compute_para_perp(module.x_err, module.y_err, module.az_tru)[1],
                "obj_tru": module.obj_tru,
                "obj_est": module.obj_est,
                "obj_err": module.obj_err
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

# Overviews

fig, axs = plt.subplots(4, 3, figsize=(18, 20))

for i, dataset in enumerate(datasets):

    data1 = dataset["para_err"]
    mean1 = np.mean(data1)
    std1 = np.std(data1)
    data2 = dataset["perp_err"]
    mean2 = np.mean(data1)
    std2 = np.std(data1)
    axs[i, 0].scatter(data1, data2, s=2, color=dgreen)
    axs[i, 0].axvline(0, color='lightgray', linewidth=1)
    axs[i, 0].axhline(0, color='lightgray', linewidth=1)
    axs[i, 0].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[i, 0].axhline(mean2, color=dred, linewidth=1)
    axs[i, 0].add_patch(
        plt.Rectangle(
            (mean1 - std1, mean2 - std2),  # Bottom-left corner
            2 * std1, # Width
            2 * std2, # Height
            color=dred,
            fill=False,
            linestyle='dashed',
            linewidth=1,
            label=f'σ = ({std1:.2f}, {std2:.2f})'
        )
    )

    axs[i, 0].set_xlim(-50, 50)
    axs[i, 0].set_ylim(-50, 50)
    axs[i, 0].set_xlabel('Δ$\parallel$')
    axs[i, 0].set_ylabel('Δ$\perp$')
    axs[i, 0].set_title(model_names[i]+"\nLocalisation residuals, 100nm")
    axs[i, 0].legend(loc='upper right')


    # Plots - zoom 40nm

    c = axs[i, 1].scatter(data1, data2, s=2, c=dataset["az_tru"], cmap='hsv')
    cbar = fig.colorbar(c)
    axs[i, 1].axvline(0, color='lightgray', linewidth=1)
    axs[i, 1].axhline(0, color='lightgray', linewidth=1)
    axs[i, 1].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[i, 1].axhline(mean2, color=dred, linewidth=1)
    axs[i, 1].add_patch(
        plt.Rectangle(
            (mean1 - std1, mean2 - std2),  # Bottom-left corner
            2 * std1, # Width
            2 * std2, # Height
            color=dred,
            fill=False,
            linestyle='dashed',
            linewidth=1,
            label=f'σ = ({std1:.2f}, {std2:.2f})'
        )
    )
    axs[i, 1].set_xlim(-5, 5)
    axs[i, 1].set_ylim(-5, 5)
    axs[i, 1].set_xlabel('Δ$\parallel$')
    axs[i, 1].set_ylabel('Δ$\perp$')
    axs[i, 1].set_title(model_names[i]+"\nLocalisation residuals, 10nm")
    axs[i, 1].legend(loc='upper right')
   
    # Plots - zoom 2σ
    c = axs[i, 2].scatter(data1, data2, s=2, c=dataset["az_tru"], cmap='hsv')
    cbar = fig.colorbar(c)
    axs[i, 2].axvline(0, color='lightgray', linewidth=1)
    axs[i, 2].axhline(0, color='lightgray', linewidth=1)
    axs[i, 2].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[i, 2].axhline(mean2, color=dred, linewidth=1)
    axs[i, 2].add_patch(
        plt.Rectangle(
            (mean1 - std1, mean2 - std2),  # Bottom-left corner
            2 * std1, # Width
            2 * std2, # Height
            color=dred,
            fill=False,
            linestyle='dashed',
            linewidth=1,
            label=f'σ = ({std1:.2f}, {std2:.2f})'
        )
    )
    axs[i, 2].set_xlim(-0.5, 0.5)
    axs[i, 2].set_ylim(-0.5, 0.5)
    axs[i, 2].set_xlabel('Δ$\parallel$')
    axs[i, 2].set_ylabel('Δ$\perp$')
    axs[i, 2].set_title(model_names[i]+"\nLocalisation residuals, 1nm")
    axs[i, 2].legend(loc='upper right')
  
  
plt.tight_layout()
#plt.show()
output_filename = "overview.png"
plt.savefig(f"{output_dir}{output_filename}", dpi=300)
plt.close()  






# Azimuth plots

fig, axs = plt.subplots(4, 5, figsize=(20, 12))

for i, dataset in enumerate(datasets):

    data = dataset["x_err"]
    mean = np.mean(data)
    std = np.std(data)
    axs[i, 0].scatter(dataset["az_tru"], data, s=2)
    axs[i, 0].axhline(0, color='lightgray', linewidth=1)
    axs[i, 0].axhline(mean, color=dred, linewidth=1, label=f'μ = {mean:.2f}\nσ = {std:.2f}')
    axs[i, 0].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
    axs[i, 0].set_ylim(-50, 50)
    axs[i, 0].set_xlabel('$\phi$, °')
    axs[i, 0].set_ylabel('Δ$\parallel$, nm')
    axs[i, 0].set_title(model_names[i]+"\nParallel localisation residuals")
    axs[i, 0].legend(loc='upper right')

    data = dataset["y_err"]
    mean = np.mean(data)
    std = np.std(data)
    axs[i, 1].scatter(dataset["az_tru"], data, s=2)
    axs[i, 1].axhline(0, color='lightgray', linewidth=1)
    axs[i, 1].axhline(mean, color=dred, linewidth=1, label=f'μ = {mean:.2f}')
    axs[i, 1].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
    axs[i, 1].set_ylim(-50, 50)
    axs[i, 1].set_xlabel('$\phi$, °')
    axs[i, 1].set_ylabel('Δ$\perp$, nm')
    axs[i, 1].set_title(model_names[i]+"\nPerpendicular localisation residuals")
    axs[i, 1].legend(loc='upper right')

    data = dataset["inc_err"]
    mean = np.mean(data)
    std = np.std(data)
    axs[i, 2].scatter(dataset["az_tru"], data, s=2)
    axs[i, 2].axhline(0, color='lightgray', linewidth=1)
    axs[i, 2].axhline(mean, color=dred, linewidth=1, label=f'μ = {mean:.2f}')
    axs[i, 2].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
    axs[i, 2].set_ylim(0, 90)
    axs[i, 2].set_xlabel('$\phi$, °')
    axs[i, 2].set_ylabel('Δθ, °')
    axs[i, 2].set_title(model_names[i]+"\nPolar residuals")
    axs[i, 2].legend(loc='upper right')

    data = dataset["az_err"]
    mean = np.mean(data)
    std = np.std(data)
    axs[i, 3].scatter(dataset["az_tru"], data, s=2)
    axs[i, 3].axhline(0, color='lightgray', linewidth=1)
    axs[i, 3].axhline(mean, color=dred, linewidth=1, label=f'μ = {mean:.2f}')
    axs[i, 3].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
    axs[i, 3].set_ylim(0, 180)
    axs[i, 3].set_xlabel('$\phi$, °')
    axs[i, 3].set_ylabel('Δ$\phi$, °')
    axs[i, 3].set_title(model_names[i]+"\nAzimuth residuals")
    axs[i, 3].legend(loc='upper right')

    data = abs(dataset["obj_err"])
    mean = np.mean(data)
    std = np.std(data)
    axs[i, 4].scatter(dataset["az_tru"], data, s=2)
    axs[i, 4].axhline(0, color='lightgray', linewidth=1)
    axs[i, 4].axhline(mean, color=dred, linewidth=1, label=f'μ = {mean:.2f}')
    axs[i, 4].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
    axs[i, 4].set_xlabel('$\phi$, °')
    axs[i, 4].set_ylabel('Δ(obj func)')
    axs[i, 4].set_title(model_names[i]+"\nObj func residuals")
    axs[i, 4].legend(loc='upper right')

plt.tight_layout()
#plt.show()
output_filename = "detail.png"
plt.savefig(f"{output_dir}{output_filename}", dpi=300)
plt.close()  


# Histograms

fig, axs = plt.subplots(4, 5, figsize=(20, 12))

for i, dataset in enumerate(datasets):

    bin_width_angles = 5

    axs[i, 0].hist(dataset["para_err"], bins=30, edgecolor='black')
    axs[i, 0].set_xlabel('Δ$\parallel$, nm')
    axs[i, 0].set_ylabel('Counts')
    axs[i, 0].set_title(model_names[i]+"\nParallel localisation residuals")

    axs[i, 1].hist(dataset["perp_err"], bins=30, edgecolor='black')
    axs[i, 1].set_xlabel('Δ$\perp$, nm')
    axs[i, 1].set_ylabel('Counts')
    axs[i, 1].set_title(model_names[i]+"\nPerpedicular localisation residuals")

    N_bins = np.arange(min(dataset["inc_err"]), max(dataset["inc_err"]) + bin_width_angles, bin_width_angles)
    axs[i, 2].hist(dataset["inc_err"], bins=N_bins, edgecolor='black')
    axs[i, 2].set_xlabel('Δθ, °')
    axs[i, 2].set_ylabel('Counts')
    axs[i, 2].set_title(model_names[i]+"\nPolar residuals")
    axs[i, 2].set_xlim(0, 90)

    N_bins = np.arange(min(dataset["az_err"]), max(dataset["az_err"]) + bin_width_angles, bin_width_angles)
    axs[i, 3].hist(dataset["az_err"], bins=N_bins, edgecolor='black')
    axs[i, 3].set_xlabel('Δ$\phi$, °')
    axs[i, 3].set_ylabel('Counts')
    axs[i, 3].set_title(model_names[i]+"\nAzimuth residuals")
    axs[i, 3].set_xlim(0, 180)

    axs[i, 4].hist(abs(dataset["obj_err"]), bins=30, edgecolor='black')
    axs[i, 4].set_xlabel('Δ(obj func)')
    axs[i, 4].set_ylabel('Counts')
    axs[i, 4].set_title(model_names[i]+"\nObj func residuals")

plt.tight_layout()
#plt.show()
output_filename = "histograms.png"
plt.savefig(f"{output_dir}{output_filename}", dpi=300)
plt.close()

