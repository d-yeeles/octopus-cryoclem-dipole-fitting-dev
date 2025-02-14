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
    "perp_err"
])

data_hinterer = data_gaussian.copy()




# Importing all the data from each results file
subdirectory = "2spot_vary_az"
for filename in os.listdir(subdirectory):
    if filename.startswith("fitting_results_inc") and filename.endswith("_gaussian.py"):

        file_path = os.path.join(subdirectory, filename)
        
        module_name = "gaussian"
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        gaussian = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(gaussian)
        
        newdata_gaussian = pd.DataFrame({
                "x_tru": gaussian.x_tru,
                "y_tru": gaussian.y_tru,
                "inc_tru": gaussian.inc_tru,
                "az_tru": gaussian.az_tru,
                "x_est": gaussian.x_est,
                "y_est": gaussian.y_est,
                "inc_est": gaussian.inc_est,
                "az_est": gaussian.az_est,
                "x_err": gaussian.x_err,
                "y_err": gaussian.y_err,
                "inc_err": gaussian.inc_err,
                "az_err": gaussian.az_err,
                "para_err": compute_para_perp(gaussian.x_err, gaussian.y_err, gaussian.az_tru)[0],
                "perp_err": compute_para_perp(gaussian.x_err, gaussian.y_err, gaussian.az_tru)[1]
            })

        if data_gaussian.empty:
            data_gaussian = newdata_gaussian
        else:
            #data_gaussian = data_gaussian.reset_index(drop=True)
            data_gaussian = pd.concat([data_gaussian, newdata_gaussian], ignore_index=True)



    elif filename.startswith("fitting_results_inc") and filename.endswith("_hinterer.py"):

        file_path = os.path.join(subdirectory, filename)
        
        module_name = "hinterer"
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        hinterer = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(hinterer)
        
        newdata_hinterer = pd.DataFrame({
                "x_tru": hinterer.x_tru,
                "y_tru": hinterer.y_tru,
                "inc_tru": hinterer.inc_tru,
                "az_tru": hinterer.az_tru,
                "x_est": hinterer.x_est,
                "y_est": hinterer.y_est,
                "inc_est": hinterer.inc_est,
                "az_est": hinterer.az_est,
                "x_err": hinterer.x_err,
                "y_err": hinterer.y_err,
                "inc_err": hinterer.inc_err,
                "az_err": hinterer.az_err,
                "para_err": compute_para_perp(hinterer.x_err, hinterer.y_err, hinterer.az_tru)[0],
                "perp_err": compute_para_perp(hinterer.x_err, hinterer.y_err, hinterer.az_tru)[1]
            })

        if data_hinterer.empty:
            data_hinterer = newdata_hinterer
        else:
            #data_hinterer = data_hinterer.reset_index(drop=True)
            data_hinterer = pd.concat([data_hinterer, newdata_hinterer], ignore_index=True)


# Convert rad to deg
angle_columns = ["inc_tru", "az_tru", "inc_est", "az_est", "inc_err", "az_err"]
data_gaussian[angle_columns] = data_gaussian[angle_columns] * 180 / np.pi
data_hinterer[angle_columns] = data_hinterer[angle_columns] * 180 / np.pi

# Account for wrapping of inc around 180 degrees:
data_gaussian["inc_err"] = np.minimum(np.abs(data_gaussian["inc_err"]), np.abs(180 - np.abs(data_gaussian["inc_err"])))
data_hinterer["inc_err"] = np.minimum(np.abs(data_hinterer["inc_err"]), np.abs(180 - np.abs(data_hinterer["inc_err"])))

# Account for wrapping of inc around 360 degrees:
data_gaussian["az_err"] = np.minimum(np.abs(data_gaussian["az_err"]), np.abs(360 - np.abs(data_gaussian["az_err"])))
data_hinterer["az_err"] = np.minimum(np.abs(data_hinterer["az_err"]), np.abs(360 - np.abs(data_hinterer["az_err"])))







output_dir = './results_plots/'



# Overviews
for inclination in [0, 23, 45, 68, 90]:

    data_gaussian_fixed_inc = data_gaussian[abs(data_gaussian['inc_tru'] - inclination) <= 5]
    data_hinterer_fixed_inc = data_hinterer[abs(data_hinterer['inc_tru'] - inclination) <= 5]

    # Count ring vs centre (optional)
    centre_definition_gaussian = data_gaussian_fixed_inc["x_err"]**2 + data_gaussian_fixed_inc["y_err"]**2 > 60**2
    frac_outliers_gaussian = centre_definition_gaussian.sum()/len(data_gaussian_fixed_inc)
    centre_definition_hinterer = data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 > 60**2
    frac_outliers_hinterer = centre_definition_hinterer.sum()/len(data_hinterer_fixed_inc)

    # Create a single figure for each plot
    fig, axs = plt.subplots(2, 3, figsize=(18, 10))

    # Hinterer - Center Only (using dyellow colormap)
    para_err_centre = data_hinterer_fixed_inc["para_err"][data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 < 60**2]
    perp_err_centre = data_hinterer_fixed_inc["perp_err"][data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 < 60**2]

    # Plot 2D histogram
    c = axs[0, 1].hist2d(para_err_centre, perp_err_centre, bins=100, cmap=dYellows)
    cbar = plt.colorbar(c[3], ax=axs[0, 1], label='Counts')

    # Draw patch over cluster 1
    if inclination == 0:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (-10, -10),  # Bottom-left corner
                20, # Width
                20, # Height
                color=dred,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 23:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (-5, -10),  # Bottom-left corner
                15, # Width
                20, # Height
                color=dred,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 45:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (-5, -10),  # Bottom-left corner
                15, # Width
                20, # Height
                color=dred,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 68:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (-5, -10),  # Bottom-left corner
                15, # Width
                20, # Height
                color=dred,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 90:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (-10, -10),  # Bottom-left corner
                20, # Width
                20, # Height
                color=dred,
                fill=False,
                linewidth=1,
            )
        )

    # Draw patch over cluster 2
    if inclination == 0:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (10, -10),  # Bottom-left corner
                20, # Width
                20, # Height
                color=dgreen,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 23:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (20, -10),  # Bottom-left corner
                20, # Width
                20, # Height
                color=dgreen,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 45:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (10, -10),  # Bottom-left corner
                10, # Width
                20, # Height
                color=dgreen,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 68:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (20, -20),  # Bottom-left corner
                17.5, # Width
                40, # Height
                color=dgreen,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 90:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (10, -10),  # Bottom-left corner
                20, # Width
                20, # Height
                color=dgreen,
                fill=False,
                linewidth=1,
            )
        )

    # Draw patch over cluster 3
    if inclination == 0:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (-30, -10),  # Bottom-left corner
                20, # Width
                20, # Height
                color=dblue,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 23:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (-40, -40),  # Bottom-left corner
                30, # Width
                80, # Height
                color=dblue,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 45:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (-35, -40),  # Bottom-left corner
                30, # Width
                80, # Height
                color=dblue,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 68:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (-25, -20),  # Bottom-left corner
                20, # Width
                40, # Height
                color=dblue,
                fill=False,
                linewidth=1,
            )
        )
    if inclination == 90:
        axs[0, 1].add_patch(
            plt.Rectangle(
                (-30, -10),  # Bottom-left corner
                20, # Width
                20, # Height
                color=dblue,
                fill=False,
                linewidth=1,
            )
        )

    # Set labels and title for the first plot
    axs[0, 1].set_xlim(-50, 50)
    axs[0, 1].set_ylim(-50, 50)
    axs[0, 1].set_xlabel('Δ$\parallel$')
    axs[0, 1].set_ylabel('Δ$\perp$')
    axs[0, 1].set_title(f'Localisation residuals histogram - Hinterer center only (Inclination: {inclination}°)')
    axs[0, 1].legend(loc='upper right')

    # Filter data based on para_err between 15 and 25
    if inclination == 0:
        filtered_data_blue = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= -30) & 
            (data_hinterer_fixed_inc["para_err"] <= -10) & 
            (data_hinterer_fixed_inc["perp_err"] >= -10) & 
            (data_hinterer_fixed_inc["perp_err"] <= 10)
        ]
        filtered_data_red = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= -10) & 
            (data_hinterer_fixed_inc["para_err"] <= 10) & 
            (data_hinterer_fixed_inc["perp_err"] >= -10) & 
            (data_hinterer_fixed_inc["perp_err"] <= 10)
        ]
        filtered_data_green = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= 10) & 
            (data_hinterer_fixed_inc["para_err"] <= 30) & 
            (data_hinterer_fixed_inc["perp_err"] >= -10) & 
            (data_hinterer_fixed_inc["perp_err"] <= 10)
        ]

    elif inclination == 23:
        filtered_data_blue = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= -40) & 
            (data_hinterer_fixed_inc["para_err"] <= -10) & 
            (data_hinterer_fixed_inc["perp_err"] >= -40) & 
            (data_hinterer_fixed_inc["perp_err"] <= 40)
        ]
        filtered_data_red = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= -5) & 
            (data_hinterer_fixed_inc["para_err"] <= 10) & 
            (data_hinterer_fixed_inc["perp_err"] >= -10) & 
            (data_hinterer_fixed_inc["perp_err"] <= 10)
        ]
        filtered_data_green = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= 20) & 
            (data_hinterer_fixed_inc["para_err"] <= 40) & 
            (data_hinterer_fixed_inc["perp_err"] >= -10) & 
            (data_hinterer_fixed_inc["perp_err"] <= 10)
        ]


    elif inclination == 45:
        filtered_data_blue = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= -35) & 
            (data_hinterer_fixed_inc["para_err"] <= -5) & 
            (data_hinterer_fixed_inc["perp_err"] >= -40) & 
            (data_hinterer_fixed_inc["perp_err"] <= 40)
        ]
        filtered_data_red = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= -5) & 
            (data_hinterer_fixed_inc["para_err"] <= 10) & 
            (data_hinterer_fixed_inc["perp_err"] >= -10) & 
            (data_hinterer_fixed_inc["perp_err"] <= 10)
        ]
        filtered_data_green = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= 10) & 
            (data_hinterer_fixed_inc["para_err"] <= 20) & 
            (data_hinterer_fixed_inc["perp_err"] >= -10) & 
            (data_hinterer_fixed_inc["perp_err"] <= 10)
        ]


    elif inclination == 68:
        filtered_data_blue = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= -25) & 
            (data_hinterer_fixed_inc["para_err"] <= -5) & 
            (data_hinterer_fixed_inc["perp_err"] >= -20) & 
            (data_hinterer_fixed_inc["perp_err"] <= 20)
        ]
        filtered_data_red = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= -5) & 
            (data_hinterer_fixed_inc["para_err"] <= 15) & 
            (data_hinterer_fixed_inc["perp_err"] >= -10) & 
            (data_hinterer_fixed_inc["perp_err"] <= 10)
        ]
        filtered_data_green = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= 20) & 
            (data_hinterer_fixed_inc["para_err"] <= 40) & 
            (data_hinterer_fixed_inc["perp_err"] >= -20) & 
            (data_hinterer_fixed_inc["perp_err"] <= 20)
        ]


    elif inclination == 90:
        filtered_data_blue = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= -30) & 
            (data_hinterer_fixed_inc["para_err"] <= -10) & 
            (data_hinterer_fixed_inc["perp_err"] >= -10) & 
            (data_hinterer_fixed_inc["perp_err"] <= 10)
        ]
        filtered_data_red = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= -10) & 
            (data_hinterer_fixed_inc["para_err"] <= 10) & 
            (data_hinterer_fixed_inc["perp_err"] >= -10) & 
            (data_hinterer_fixed_inc["perp_err"] <= 10)
        ]
        filtered_data_green = data_hinterer_fixed_inc[
            (data_hinterer_fixed_inc["para_err"] >= 10) & 
            (data_hinterer_fixed_inc["para_err"] <= 30) & 
            (data_hinterer_fixed_inc["perp_err"] >= -10) & 
            (data_hinterer_fixed_inc["perp_err"] <= 10)
        ]


    mean_blue = np.mean(filtered_data_blue["az_tru"])    
    mean_red = np.mean(filtered_data_red["az_tru"])    
    mean_green = np.mean(filtered_data_green["az_tru"])    
    std_blue = np.std(filtered_data_blue["az_tru"])    
    std_red = np.std(filtered_data_red["az_tru"])    
    std_green = np.std(filtered_data_green["az_tru"])    

    # Second plot: Histogram of 'az_tru'
    axs[1, 0].hist(filtered_data_blue["az_tru"], bins=36, color=dblue, edgecolor='black')
    axs[1, 0].axvline(mean_blue, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean_blue:.2f}')
    axs[1, 0].fill_betweenx(
        y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
        x1=mean_blue - std_blue,
        x2=mean_blue + std_blue,
        color=dblue,
        alpha=0.1,
        label=f'σ = {std_blue:.2f}'
    )
    axs[1, 0].set_title('Distribution of ϕ in blue cluster')
    axs[1, 0].set_xlabel('az (deg)')
    axs[1, 0].set_ylabel('Frequency')

    # Second plot: Histogram of 'az_tru'
    axs[1, 1].hist(filtered_data_red["az_tru"], bins=36, color=dred, edgecolor='black')
    axs[1, 1].axvline(mean_red, color=dred, linestyle='dashed', linewidth=1, label=f'μ = {mean_red:.2f}')
    axs[1, 1].fill_betweenx(
        y=[0, axs[1, 1].get_ylim()[1]],  # Cover full y-axis range
        x1=mean_red - std_red,
        x2=mean_red + std_red,
        color=dred,
        alpha=0.1,
        label=f'σ = {std_red:.2f}'
    )
    axs[1, 1].set_title('Distribution of ϕ in red cluster')
    axs[1, 1].set_xlabel('az (deg)')
    axs[1, 1].set_ylabel('Frequency')

    # Second plot: Histogram of 'az_tru'
    axs[1, 2].hist(filtered_data_green["az_tru"], bins=36, color=dgreen, edgecolor='black')
    axs[1, 2].axvline(mean_green, color=dgreen, linestyle='dashed', linewidth=1, label=f'μ = {mean_red:.2f}')
    axs[1, 2].fill_betweenx(
        y=[0, axs[1, 2].get_ylim()[1]],  # Cover full y-axis range
        x1=mean_green - std_green,
        x2=mean_green + std_green,
        color=dgreen,
        alpha=0.1,
        label=f'σ = {std_green:.2f}'
    )
    axs[1, 2].set_title('Distribution of ϕ in green cluster')
    axs[1, 2].set_xlabel('az (deg)')
    axs[1, 2].set_ylabel('Frequency')





    plt.tight_layout()
#    if inclination == 68i:
    output_filename = f"clusters_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}")
#    else:
#        plt.show()
    plt.close()
