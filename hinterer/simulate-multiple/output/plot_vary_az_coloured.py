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

# Create custom colormaps for hists
dBlues = LinearSegmentedColormap.from_list('dblue_to_white', [(1, 1, 1), dblue], N=100)
dYellows = LinearSegmentedColormap.from_list('dyellow_to_white', [(1, 1, 1), dyellow], N=100)


# Initialise empty dataframes
data_hinterer1 = pd.DataFrame(columns=[
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
    "perp_err"
])

data_hinterer2 = data_hinterer1.copy()
data_hinterer3 = data_hinterer2.copy()




# Importing all the data from each results file
subdirectory = "2spot_vary_az"
for filename in os.listdir(subdirectory):
    if filename.startswith("fitting_results_inc") and filename.endswith("_hinterer.py"):
        file_path = os.path.join(subdirectory, filename)
        
        module_name = "hinterer1"
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        hinterer1 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(hinterer1)
        
        newdata_hinterer1 = pd.DataFrame({
                "x_tru": hinterer1.x_tru,
                "y_tru": hinterer1.y_tru,
                "inc_tru": hinterer1.inc_tru,
                "az_tru": hinterer1.az_tru,
                "x_est": hinterer1.x_est,
                "y_est": hinterer1.y_est,
                "inc_est": hinterer1.inc_est,
                "az_est": hinterer1.az_est,
                "x_err": hinterer1.x_err,
                "y_err": hinterer1.y_err,
                "inc_err": hinterer1.inc_err,
                "az_err": hinterer1.az_err,
                "para_err": compute_para_perp(hinterer1.x_err, hinterer1.y_err, hinterer1.az_tru)[0],
                "perp_err": compute_para_perp(hinterer1.x_err, hinterer1.y_err, hinterer1.az_tru)[1]
            })

        if data_hinterer1.empty:
            data_hinterer1 = newdata_hinterer1
        else:
            data_hinterer1 = pd.concat([data_hinterer1, newdata_hinterer1], ignore_index=True)



# Importing all the data from each results file
subdirectory = "2spot_inc68_v2"
for filename in os.listdir(subdirectory):
    if filename.startswith("redo1_fitting_results_inc"):

        file_path = os.path.join(subdirectory, filename)

        module_name = "hinterer2"
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        hinterer2 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(hinterer2)

        newdata_hinterer2 = pd.DataFrame({
                "x_tru": hinterer2.x_tru,
                "y_tru": hinterer2.y_tru,
                "inc_tru": hinterer2.inc_tru,
                "az_tru": hinterer2.az_tru,
                "x_est": hinterer2.x_est,
                "y_est": hinterer2.y_est,
                "inc_est": hinterer2.inc_est,
                "az_est": hinterer2.az_est,
                "x_err": hinterer2.x_err,
                "y_err": hinterer2.y_err,
                "inc_err": hinterer2.inc_err,
                "az_err": hinterer2.az_err,
                "para_err": compute_para_perp(hinterer2.x_err, hinterer2.y_err, hinterer2.az_tru)[0],
                "perp_err": compute_para_perp(hinterer2.x_err, hinterer2.y_err, hinterer2.az_tru)[1]
            })

        if data_hinterer2.empty:
            data_hinterer2 = newdata_hinterer2
        else:
            data_hinterer2 = pd.concat([data_hinterer2, newdata_hinterer2], ignore_index=True)


subdirectory = "2spot_inc68_v2"
for filename in os.listdir(subdirectory):
    if filename.startswith("redo2_fitting_results_inc"):

        file_path = os.path.join(subdirectory, filename)

        module_name = "hinterer3"
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        hinterer3 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(hinterer3)

        newdata_hinterer3 = pd.DataFrame({
                "x_tru": hinterer3.x_tru,
                "y_tru": hinterer3.y_tru,
                "inc_tru": hinterer3.inc_tru,
                "az_tru": hinterer3.az_tru,
                "x_est": hinterer3.x_est,
                "y_est": hinterer3.y_est,
                "inc_est": hinterer3.inc_est,
                "az_est": hinterer3.az_est,
                "x_err": hinterer3.x_err,
                "y_err": hinterer3.y_err,
                "inc_err": hinterer3.inc_err,
                "az_err": hinterer3.az_err,
                "para_err": compute_para_perp(hinterer3.x_err, hinterer3.y_err, hinterer3.az_tru)[0],
                "perp_err": compute_para_perp(hinterer3.x_err, hinterer3.y_err, hinterer3.az_tru)[1]
            })

        if data_hinterer3.empty:
            data_hinterer3 = newdata_hinterer3
        else:
            data_hinterer3 = pd.concat([data_hinterer3, newdata_hinterer3], ignore_index=True)





# Convert rad to deg
angle_columns = ["inc_tru", "az_tru", "inc_est", "az_est", "inc_err", "az_err"]
data_hinterer1[angle_columns] = data_hinterer1[angle_columns] * 180 / np.pi
data_hinterer2[angle_columns] = data_hinterer2[angle_columns] * 180 / np.pi
data_hinterer3[angle_columns] = data_hinterer3[angle_columns] * 180 / np.pi

# Account for wrapping of inc around 180 degrees:
data_hinterer1["inc_err"] = np.minimum(np.abs(data_hinterer1["inc_err"]), np.abs(180 - np.abs(data_hinterer1["inc_err"])))
data_hinterer2["inc_err"] = np.minimum(np.abs(data_hinterer2["inc_err"]), np.abs(180 - np.abs(data_hinterer2["inc_err"])))
data_hinterer3["inc_err"] = np.minimum(np.abs(data_hinterer3["inc_err"]), np.abs(180 - np.abs(data_hinterer3["inc_err"])))

# Account for wrapping of inc around 360 degrees:
data_hinterer1["az_err"] = np.minimum(np.abs(data_hinterer1["az_err"]), np.abs(360 - np.abs(data_hinterer1["az_err"])))
data_hinterer2["az_err"] = np.minimum(np.abs(data_hinterer2["az_err"]), np.abs(360 - np.abs(data_hinterer2["az_err"])))
data_hinterer3["az_err"] = np.minimum(np.abs(data_hinterer3["az_err"]), np.abs(360 - np.abs(data_hinterer3["az_err"])))






output_dir = './results_plots/'



# Overviews
fig, axs = plt.subplots(1, 3, figsize=(36, 6))

for inc_index, inclination in enumerate([68]):

    data_hinterer1_fixed_inc = data_hinterer1[abs(data_hinterer1['inc_tru'] - inclination) <= 5]
    data_hinterer2_fixed_inc = data_hinterer2[abs(data_hinterer2['inc_tru'] - inclination) <= 5]
    data_hinterer3_fixed_inc = data_hinterer3[abs(data_hinterer3['inc_tru'] - inclination) <= 5]


    # Count ring vs centre (optional)
    centre_definition_hinterer1 = data_hinterer1_fixed_inc["x_err"]**2 + data_hinterer1_fixed_inc["y_err"]**2 > 60**2
    frac_outliers_hinterer1 = centre_definition_hinterer1.sum()/len(data_hinterer1_fixed_inc)
    centre_definition_hinterer2 = data_hinterer2_fixed_inc["x_err"]**2 + data_hinterer2_fixed_inc["y_err"]**2 > 60**2
    frac_outliers_hinterer2 = centre_definition_hinterer2.sum()/len(data_hinterer2_fixed_inc)
    centre_definition_hinterer3 = data_hinterer3_fixed_inc["x_err"]**2 + data_hinterer3_fixed_inc["y_err"]**2 > 60**2
    frac_outliers_hinterer3 = centre_definition_hinterer3.sum()/len(data_hinterer3_fixed_inc)

    # Hinterer old (centre only)
    para_err_centre1 = data_hinterer1_fixed_inc["para_err"][data_hinterer1_fixed_inc["x_err"]**2 + data_hinterer1_fixed_inc["y_err"]**2 < 60**2]
    perp_err_centre1 = data_hinterer1_fixed_inc["perp_err"][data_hinterer1_fixed_inc["x_err"]**2 + data_hinterer1_fixed_inc["y_err"]**2 < 60**2]
    az_tru_centre1 = data_hinterer1_fixed_inc["az_tru"][data_hinterer1_fixed_inc["x_err"]**2 + data_hinterer1_fixed_inc["y_err"]**2 < 60**2]
    # Hinterer new (centre only)
    para_err_centre2 = data_hinterer2_fixed_inc["para_err"][data_hinterer2_fixed_inc["x_err"]**2 + data_hinterer2_fixed_inc["y_err"]**2 < 60**2]
    perp_err_centre2 = data_hinterer2_fixed_inc["perp_err"][data_hinterer2_fixed_inc["x_err"]**2 + data_hinterer2_fixed_inc["y_err"]**2 < 60**2]
    az_tru_centre2 = data_hinterer2_fixed_inc["az_tru"][data_hinterer2_fixed_inc["x_err"]**2 + data_hinterer2_fixed_inc["y_err"]**2 < 60**2]
        # Hinterer new (centre only)
    para_err_centre3 = data_hinterer3_fixed_inc["para_err"][data_hinterer3_fixed_inc["x_err"]**2 + data_hinterer3_fixed_inc["y_err"]**2 < 60**2]
    perp_err_centre3 = data_hinterer3_fixed_inc["perp_err"][data_hinterer3_fixed_inc["x_err"]**2 + data_hinterer3_fixed_inc["y_err"]**2 < 60**2]
    az_tru_centre3 = data_hinterer3_fixed_inc["az_tru"][data_hinterer3_fixed_inc["x_err"]**2 + data_hinterer3_fixed_inc["y_err"]**2 < 60**2]

    # Inc = 0
#    c = axs[inc_index].hist2d(para_err_centre, perp_err_centre, bins=100, cmap=dYellows)
#    cbar = fig.colorbar(c[3], ax=axs[inc_index], label='Counts')
    scatter = axs[0].scatter(para_err_centre1, perp_err_centre1, c=az_tru_centre1, cmap='cool', s=1)#, alpha=0.1)
    cbar = fig.colorbar(scatter, ax=axs[0], label='ϕ')
    axs[0].set_xlim(-50, 50)
    axs[0].set_ylim(-50, 50)
    axs[0].set_xlabel('Δ$\parallel$')
    axs[0].set_ylabel('Δ$\perp$')
    axs[0].set_title(f'Old version (50 per ϕ)\nθ = {inclination}°')
    axs[0].set_aspect('equal')

    scatter = axs[1].scatter(para_err_centre2, perp_err_centre2, c=az_tru_centre2, cmap='cool', s=1)#, alpha=0.1)
    cbar = fig.colorbar(scatter, ax=axs[1], label='ϕ')
    axs[1].set_xlim(-50, 50)
    axs[1].set_ylim(-50, 50)
    axs[1].set_xlabel('Δ$\parallel$')
    axs[1].set_ylabel('Δ$\perp$')
    axs[1].set_title(f'New version (8 per ϕ)\nθ = {inclination}°')
    axs[1].set_aspect('equal')

    scatter = axs[2].scatter(para_err_centre3, perp_err_centre3, c=az_tru_centre3, cmap='cool', s=1)#, alpha=0.1)
    cbar = fig.colorbar(scatter, ax=axs[2], label='ϕ')
    axs[2].set_xlim(-50, 50)
    axs[2].set_ylim(-50, 50)
    axs[2].set_xlabel('Δ$\parallel$')
    axs[2].set_ylabel('Δ$\perp$')
    axs[2].set_title(f'New version (8 per ϕ)\nθ = {inclination}°')
    axs[2].set_aspect('equal')

plt.tight_layout()
output_filename = f"coloured_azimuths_inc{round(inclination)}_68only.png"
plt.savefig(f"{output_dir}{output_filename}")
plt.show()
plt.close()
