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
subdirectory = "../2spot_vary_az"
for filename in os.listdir(subdirectory):
#    if filename.startswith("fitting_results_inc") and filename.endswith("_gaussian.py"):
    if filename.startswith("fitting_results_gaussian_redo"):

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



#    elif filename.startswith("fitting_results_inc") and filename.endswith("_hinterer.py"):
subdirectory = "./"
for filename in os.listdir(subdirectory):
    if filename.startswith("fitting_results_swarm_inc0_4"):

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
                "perp_err": compute_para_perp(hinterer.x_err, hinterer.y_err, hinterer.az_tru)[1],
                "obj_tru": hinterer.obj_tru,
                "obj_est": hinterer.obj_est,
                "obj_err": hinterer.obj_err

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
#data_gaussian["inc_err"] = np.minimum(np.abs(data_gaussian["inc_err"]), np.abs(180 - np.abs(data_gaussian["inc_err"])))
#data_hinterer["inc_err"] = np.minimum(np.abs(data_hinterer["inc_err"]), np.abs(180 - np.abs(data_hinterer["inc_err"])))
data_gaussian["inc_err"] = np.mod(data_gaussian["inc_err"], 180)
data_hinterer["inc_err"] = np.mod(data_hinterer["inc_err"], 180)
data_gaussian["inc_err"] = np.minimum(np.abs(data_gaussian["inc_err"]), 180 - np.abs(data_gaussian["inc_err"]))
data_hinterer["inc_err"] = np.minimum(np.abs(data_hinterer["inc_err"]), 180 - np.abs(data_hinterer["inc_err"]))

# Account for wrapping of az around 360 degrees:
#data_gaussian["az_err"] = np.minimum(np.abs(data_gaussian["az_err"]), np.abs(360 - np.abs(data_gaussian["az_err"])))
#data_hinterer["az_err"] = np.minimum(np.abs(data_hinterer["az_err"]), 360 - np.abs(data_hinterer["az_err"]))
data_gaussian["az_err"] = np.mod(data_gaussian["az_err"], 360)
data_hinterer["az_err"] = np.mod(data_hinterer["az_err"], 360)
data_gaussian["az_err"] = np.minimum(np.abs(data_gaussian["az_err"]), 360 - np.abs(data_gaussian["az_err"]))
data_hinterer["az_err"] = np.minimum(np.abs(data_hinterer["az_err"]), 360 - np.abs(data_hinterer["az_err"]))





output_dir = './'

# Overviews

# Generate plots for each fixed inclination
for inclination in [0, 23, 45, 68, 90]:

    #data_gaussian_fixed_inc = data_gaussian[round(data_gaussian['inc_tru'] / 5) * 5 == inclination]
    #data_hinterer_fixed_inc = data_hinterer[round(data_hinterer['inc_tru'] / 5) * 5 == inclination]
    data_gaussian_fixed_inc = data_gaussian[abs(data_gaussian['inc_tru'] - inclination) <= 5]
    data_hinterer_fixed_inc = data_hinterer[abs(data_hinterer['inc_tru'] - inclination) <= 5]

    # Count ring vs centre
    centre_definition_gaussian = data_gaussian_fixed_inc["x_err"]**2 + data_gaussian_fixed_inc["y_err"]**2 > 60**2
    frac_outliers_gaussian = centre_definition_gaussian.sum()/len(data_gaussian_fixed_inc)
    centre_definition_hinterer = data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 > 60**2
    frac_outliers_hinterer = centre_definition_hinterer.sum()/len(data_hinterer_fixed_inc)

    # Load and display background image
    background_image = plt.imread(f'../animation_lowres/vary_inc/sim_inc{round(inclination):03d}_az0_run1.png')

    # Gaussian
    # Plots - full
    fig, axs = plt.subplots(2, 3, figsize=(18, 10))
    
    para_err_full = data_gaussian_fixed_inc["para_err"]
    perp_err_full = data_gaussian_fixed_inc["perp_err"]
    para_err_centre = data_gaussian_fixed_inc["para_err"][data_gaussian_fixed_inc["x_err"]**2 + data_gaussian_fixed_inc["y_err"]**2 < 60**2]
    perp_err_centre = data_gaussian_fixed_inc["perp_err"][data_gaussian_fixed_inc["x_err"]**2 + data_gaussian_fixed_inc["y_err"]**2 < 60**2]
    para_err_ring = data_gaussian_fixed_inc["para_err"][data_gaussian_fixed_inc["x_err"]**2 + data_gaussian_fixed_inc["y_err"]**2 > 60**2]
    perp_err_ring = data_gaussian_fixed_inc["perp_err"][data_gaussian_fixed_inc["x_err"]**2 + data_gaussian_fixed_inc["y_err"]**2 > 60**2]
    
    mean1 = np.mean(para_err_full)
    std1 = np.std(para_err_full)
    mean2 = np.mean(perp_err_full)
    std2 = np.std(perp_err_full)
    axs[0, 0].scatter(para_err_full, perp_err_full, s=2, color=dblue)
    axs[0, 0].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[0, 0].axhline(mean2, color=dred, linewidth=1)
    axs[0, 0].add_patch(
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

    axs[0, 0].set_facecolor('black')
    axs[0, 0].imshow(background_image, cmap='gray', extent=(-500, 500, -500, 500), aspect='auto', zorder=-1)

    axs[0, 0].set_xlim(-700, 700)
    axs[0, 0].set_ylim(-700, 700)
    axs[0, 0].set_xlabel('Δ$\parallel$')
    axs[0, 0].set_ylabel('Δ$\perp$')
    axs[0, 0].set_title("Overview of localisation residuals")
    axs[0, 0].set_title(f"Overview of localisation residuals\n{np.round(frac_outliers_gaussian*100, 2)}% in ring")
    axs[0, 0].legend(loc='upper right')


    # Plots - centre only

    mean1 = np.mean(para_err_centre)
    std1 = np.std(para_err_centre)
    mean2 = np.mean(perp_err_centre)
    std2 = np.std(perp_err_centre)
#    c = axs[0, 1].hist2d(para_err_centre, perp_err_centre, bins=100, cmap=dBlues)
#    cbar = fig.colorbar(c[3], ax=axs[0, 1], label='Counts')
    axs[0, 1].scatter(para_err_centre, perp_err_centre, s=2, c=dblue)
    axs[0, 1].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[0, 1].axhline(mean2, color=dred, linewidth=1)
    axs[0, 1].add_patch(
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
    axs[0, 1].set_xlim(-50, 50)
    axs[0, 1].set_ylim(-50, 50)
    axs[0, 1].set_xlabel('Δ$\parallel$')
    axs[0, 1].set_ylabel('Δ$\perp$')
    axs[0, 1].set_title('Localisation residuals histogram - centre only')
    axs[0, 1].legend(loc='upper right')
   
    # Plots - ring only
    mean1 = np.mean(para_err_ring)
    std1 = np.std(para_err_ring)
    mean2 = np.mean(perp_err_ring)
    std2 = np.std(perp_err_ring)
    c = axs[0, 2].hist2d(para_err_ring, perp_err_ring, bins=100, cmap=dBlues)
    cbar = fig.colorbar(c[3], ax=axs[0, 2], label='Counts')
    cbar.set_ticks(np.arange(int(np.min(c[0])), int(np.max(c[0])) + 1, 1))  # Only show integers
    axs[0, 2].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[0, 2].axhline(mean2, color=dred, linewidth=1)
    axs[0, 2].add_patch(
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
    axs[0, 2].set_xlim(-700, 700)
    axs[0, 2].set_ylim(-700, 700)
    axs[0, 2].set_xlabel('Δ$\parallel$')
    axs[0, 2].set_ylabel('Δ$\perp$')
    axs[0, 2].set_title('Localisation residuals histogram - ring only')
    axs[0, 2].legend(loc='upper right')

 
    # Hinterer
    
    para_err_full = data_hinterer_fixed_inc["para_err"]
    perp_err_full = data_hinterer_fixed_inc["perp_err"]
    para_err_centre = data_hinterer_fixed_inc["para_err"][data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 < 60**2]
    perp_err_centre = data_hinterer_fixed_inc["perp_err"][data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 < 60**2]
    para_err_ring = data_hinterer_fixed_inc["para_err"][data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 > 60**2]
    perp_err_ring = data_hinterer_fixed_inc["perp_err"][data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 > 60**2]
    
    # Plots - full
    mean1 = np.mean(para_err_full)
    std1 = np.std(para_err_full)
    mean2 = np.mean(perp_err_full)
    std2 = np.std(perp_err_full)
    axs[1, 0].scatter(para_err_full, perp_err_full, s=2, color=dyellow)
    axs[1, 0].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[1, 0].axhline(mean2, color=dred, linewidth=1)
    axs[1, 0].add_patch(
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
    axs[1, 0].set_facecolor('black')
    axs[1, 0].imshow(background_image, cmap='gray', extent=(-500, 500, -500, 500), aspect='auto', zorder=-1)

    axs[1, 0].set_xlim(-700, 700)
    axs[1, 0].set_ylim(-700, 700)
    axs[1, 0].set_xlabel('Δ$\parallel$')
    axs[1, 0].set_ylabel('Δ$\perp$')
    axs[1, 0].set_title("Overview of localisation residuals")
    axs[1, 0].set_title(f"Overview of localisation residuals\n{np.round(frac_outliers_hinterer*100, 2)}% in ring")
    axs[1, 0].legend(loc='upper right')

    # Plots - centre only
    
    obj_func_err = (data_hinterer_fixed_inc["obj_tru"][data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 < 60**2])

    mean1 = np.mean(para_err_centre)
    std1 = np.std(para_err_centre)
    mean2 = np.mean(perp_err_centre)
    std2 = np.std(perp_err_centre)
    #c = axs[1, 1].hist2d(para_err_centre, perp_err_centre, bins=100, cmap=dYellows)
    #cbar = fig.colorbar(c[3], ax=axs[1, 1], label='Counts')
    #axs[1, 1].scatter(para_err_centre, perp_err_centre, s=2, color=dyellow, alpha=0.75)
    c = axs[1, 1].scatter(para_err_centre, perp_err_centre, s=2, c=obj_func_err, cmap='plasma')
    cbar = fig.colorbar(c)
    axs[1, 1].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[1, 1].axhline(mean2, color=dred, linewidth=1)
    axs[1, 1].add_patch(
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
    axs[1, 1].set_xlim(-50, 50)
    axs[1, 1].set_ylim(-50, 50)
    axs[1, 1].set_xlabel('Δ$\parallel$')
    axs[1, 1].set_ylabel('Δ$\perp$')
    axs[1, 1].set_title('Localisation residuals histogram - centre only')
    axs[1, 1].legend(loc='upper right')

    # Plots - ring only
    mean1 = np.mean(para_err_ring)
    std1 = np.std(para_err_ring)
    mean2 = np.mean(perp_err_ring)
    std2 = np.std(perp_err_ring)
    c = axs[1, 2].hist2d(para_err_ring, perp_err_ring, bins=100, cmap=dYellows)
    cbar = fig.colorbar(c[3], ax=axs[1, 2], label='Counts')
    cbar.set_ticks(np.arange(int(np.min(c[0])), int(np.max(c[0])) + 1, 1))  # Only show integers
    axs[1, 2].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[1, 2].axhline(mean2, color=dred, linewidth=1)
    axs[1, 2].add_patch(
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
    axs[1, 2].set_xlim(-700, 700)
    axs[1, 2].set_ylim(-700, 700)
    axs[1, 2].set_xlabel('Δ$\parallel$')
    axs[1, 2].set_ylabel('Δ$\perp$')
    axs[1, 2].set_title('Localisation residuals histogram - ring only')
    axs[1, 2].legend(loc='upper right')

    plt.tight_layout()
    #plt.show()
    output_filename = f"overview_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}")
    plt.close()     
     







# Histograms etc.

# Generate plots for each fixed inclination
for inclination in [0, 23, 45, 68, 90]:

    for index, data_model in enumerate([data_gaussian, data_hinterer]):

        plot_colour = [dblue, dyellow][index]
        hist_colour = [dBlues, dYellows][index]
        plot_label = ['gaussian', 'hinterer'][index]

        data_model_fixed_inc = data_model[abs(data_model['inc_tru'] - inclination) <= 5]




        # Gaussian histograms
        
        # Plots - inc/az gaussian
        fig, axs = plt.subplots(2, 3, figsize=(18, 10))
        
        # Scatter plots
        az_tru = data_model_fixed_inc["az_tru"][data_model_fixed_inc["x_err"]**2 + data_model_fixed_inc["y_err"]**2 < 60**2]
        para_err = data_model_fixed_inc["para_err"][data_model_fixed_inc["x_err"]**2 + data_model_fixed_inc["y_err"]**2 < 60**2]
        perp_err = data_model_fixed_inc["perp_err"][data_model_fixed_inc["x_err"]**2 + data_model_fixed_inc["y_err"]**2 < 60**2]
        az_err = data_model_fixed_inc["az_err"][data_model_fixed_inc["x_err"]**2 + data_model_fixed_inc["y_err"]**2 < 60**2]
        
        axs[0, 0].scatter(az_tru, para_err, s=10, alpha=0.1, color=plot_colour, label=plot_label)
        axs[0, 0].axhline(0, color=dred, linewidth=2)
        axs[0, 0].set_ylim(-50, 50)
        axs[0, 0].set_xlabel('φ, °')
        axs[0, 0].set_ylabel('Δ$\parallel$, nm')
        axs[0, 0].set_title("Localisation residuals, parallel - centre only")
        axs[0, 0].legend(loc='upper right')
        
        axs[0, 1].scatter(az_tru, perp_err, s=10, alpha=0.1, color=plot_colour, label=plot_label)
        axs[0, 1].axhline(0, color=dred, linewidth=2)
        axs[0, 1].set_ylim(-50, 50)
        axs[0, 1].set_xlabel('φ, °')
        axs[0, 1].set_ylabel('Δ$\perp$, nm')
        axs[0, 1].set_title("Localisation residuals, perpendicular - centre only")
        axs[0, 1].legend(loc='upper right')
        
        axs[0, 2].scatter(az_tru, az_err, s=10, alpha=0.1, color=plot_colour, label=plot_label)
        axs[0, 2].axhline(0, color=dred, linewidth=2)
        axs[0, 2].set_xlabel('φ, °')
        axs[0, 2].set_ylabel('Δφ, °')
        axs[0, 2].set_title("Localisation residuals, azimuth - centre only")
        axs[0, 2].legend(loc='upper right')
        
        
        # Histogram of ∥ errors
        mean = np.mean(para_err)
        std = np.std(para_err)
        axs[1, 0].hist(para_err, bins=20, alpha=0.7, color=plot_colour, label=plot_label)
        axs[1, 0].axvline(mean, color=plot_colour, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
        axs[1, 0].fill_betweenx(
            y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
            x1=mean - std,
            x2=mean + std,
            color=plot_colour,
            alpha=0.1,
            label=f'σ = {std:.2f}'
        )
        axs[1, 0].set_xlim(-50, 50)
        axs[1, 0].set_xlabel('Δ$\parallel$, nm')
        axs[1, 0].set_ylabel('Frequency')
        axs[1, 0].set_title("Parallel residual histogram - centre only")
        axs[1, 0].legend(loc='upper right')
        
        # Histogram of ⟂ errors
        mean = np.mean(perp_err)
        std = np.std(perp_err)
        axs[1, 1].hist(perp_err, bins=20, alpha=0.7, color=plot_colour, label=plot_label)
        axs[1, 1].axvline(mean, color=plot_colour, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
        axs[1, 1].fill_betweenx(
            y=[0, axs[1, 1].get_ylim()[1]],  # Cover full y-axis range
            x1=mean - std,
            x2=mean + std,
            color=plot_colour,
            alpha=0.1,
            label=f'σ = {std:.2f}'
        )
        axs[1, 1].set_xlim(-50, 50)
        axs[1, 1].set_xlabel('Δ$\perp$, nm')
        axs[1, 1].set_ylabel('Frequency')
        axs[1, 1].set_title("Perpendicular residual histogram - centre only")
        axs[1, 1].legend(loc='upper right')
        
        # Histogram of az errors
        mean = np.mean(az_err)
        std = np.std(az_err)
        axs[1, 2].hist(az_err, bins=20, alpha=0.7, color=plot_colour, label=plot_label)
        axs[1, 2].axvline(mean, color=plot_colour, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
        axs[1, 2].fill_betweenx(
            y=[0, axs[1, 2].get_ylim()[1]],  # Cover full y-axis range
            x1=mean - std,
            x2=mean + std,
            color=plot_colour,
            alpha=0.1,
            label=f'σ = {std:.2f}'
        )
        axs[1, 2].set_xlabel('Δφ, °')
        axs[1, 2].set_ylabel('Frequency')
        axs[1, 2].set_title("Azimuthal residual histogram, centre only")
        axs[1, 2].legend(loc='upper right')
        
        plt.tight_layout()

        #plt.show()
        output_filename = f"{plot_label}_inc{round(inclination)}.png"
        plt.savefig(f"{output_dir}{output_filename}")
        plt.close()     

