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
dpink = default_colors[6]

# Create custom colormaps for hists
dBlues = LinearSegmentedColormap.from_list('dblue_to_white', [(1, 1, 1), dblue], N=100)
dYellows = LinearSegmentedColormap.from_list('dyellow_to_white', [(1, 1, 1), dyellow], N=100)
dGreens = LinearSegmentedColormap.from_list('dgreen_to_white', [(1, 1, 1), dgreen], N=100)
dPinks = LinearSegmentedColormap.from_list('dpink_to_white', [(1, 1, 1), dpink], N=100)


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
data_givenpolar = data_gaussian.copy()
data_givenorient = data_gaussian.copy()




# Importing all the data from each results file
subdirectory = "2spot_vary_az"
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
    elif filename.startswith("fitting_results_hinterer_redo"):

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

    elif filename.startswith("fitting_results_correct_polar"):

        file_path = os.path.join(subdirectory, filename)

        module_name = "givenpolar"
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        givenpolar = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(givenpolar)

        newdata_givenpolar = pd.DataFrame({
                "x_tru": givenpolar.x_tru,
                "y_tru": givenpolar.y_tru,
                "inc_tru": givenpolar.inc_tru,
                "az_tru": givenpolar.az_tru,
                "x_est": givenpolar.x_est,
                "y_est": givenpolar.y_est,
                "inc_est": givenpolar.inc_est,
                "az_est": givenpolar.az_est,
                "x_err": givenpolar.x_err,
                "y_err": givenpolar.y_err,
                "inc_err": givenpolar.inc_err,
                "az_err": givenpolar.az_err,
                "para_err": compute_para_perp(givenpolar.x_err, givenpolar.y_err, givenpolar.az_tru)[0],
                "perp_err": compute_para_perp(givenpolar.x_err, givenpolar.y_err, givenpolar.az_tru)[1]
            })

        if data_givenpolar.empty:
            data_givenpolar = newdata_givenpolar
        else:
            #data_givenpolar = data_givenpolar.reset_index(drop=True)
            data_givenpolar = pd.concat([data_givenpolar, newdata_givenpolar], ignore_index=True)

    elif filename.startswith("fitting_results_correct_orient"):

        file_path = os.path.join(subdirectory, filename)

        module_name = "givenorient"
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        givenorient = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(givenorient)

        newdata_givenorient = pd.DataFrame({
                "x_tru": givenorient.x_tru,
                "y_tru": givenorient.y_tru,
                "inc_tru": givenorient.inc_tru,
                "az_tru": givenorient.az_tru,
                "x_est": givenorient.x_est,
                "y_est": givenorient.y_est,
                "inc_est": givenorient.inc_est,
                "az_est": givenorient.az_est,
                "x_err": givenorient.x_err,
                "y_err": givenorient.y_err,
                "inc_err": givenorient.inc_err,
                "az_err": givenorient.az_err,
                "para_err": compute_para_perp(givenorient.x_err, givenorient.y_err, givenorient.az_tru)[0],
                "perp_err": compute_para_perp(givenorient.x_err, givenorient.y_err, givenorient.az_tru)[1]
            })

        if data_givenorient.empty:
            data_givenorient = newdata_givenorient
        else:
            #data_givenorient = data_givenorient.reset_index(drop=True)
            data_givenorient = pd.concat([data_givenorient, newdata_givenorient], ignore_index=True)


# -----------------------


# Convert rad to deg
angle_columns = ["inc_tru", "az_tru", "inc_est", "az_est", "inc_err", "az_err"]
data_gaussian[angle_columns] = data_gaussian[angle_columns] * 180 / np.pi
data_hinterer[angle_columns] = data_hinterer[angle_columns] * 180 / np.pi
data_givenpolar[angle_columns] = data_givenpolar[angle_columns] * 180 / np.pi
data_givenorient[angle_columns] = data_givenorient[angle_columns] * 180 / np.pi

# Account for wrapping of inc around 180 degrees:
data_gaussian["inc_err"] = np.mod(data_gaussian["inc_err"], 180)
data_hinterer["inc_err"] = np.mod(data_hinterer["inc_err"], 180)
data_givenpolar["inc_err"] = np.mod(data_givenpolar["inc_err"], 180)
data_givenorient["inc_err"] = np.mod(data_givenorient["inc_err"], 180)
data_gaussian["inc_err"] = np.minimum(np.abs(data_gaussian["inc_err"]), 180 - np.abs(data_gaussian["inc_err"]))
data_hinterer["inc_err"] = np.minimum(np.abs(data_hinterer["inc_err"]), 180 - np.abs(data_hinterer["inc_err"]))
data_givenpolar["inc_err"] = np.minimum(np.abs(data_givenpolar["inc_err"]), 180 - np.abs(data_givenpolar["inc_err"]))
data_givenorient["inc_err"] = np.minimum(np.abs(data_givenorient["inc_err"]), 180 - np.abs(data_givenorient["inc_err"]))

# Account for wrapping of az around 360 degrees:
data_gaussian["az_err"] = np.mod(data_gaussian["az_err"], 360)
data_hinterer["az_err"] = np.mod(data_hinterer["az_err"], 360)
data_givenpolar["az_err"] = np.mod(data_givenpolar["az_err"], 360)
data_givenorient["az_err"] = np.mod(data_givenorient["az_err"], 360)
data_gaussian["az_err"] = np.minimum(np.abs(data_gaussian["az_err"]), 360 - np.abs(data_gaussian["az_err"]))
data_hinterer["az_err"] = np.minimum(np.abs(data_hinterer["az_err"]), 360 - np.abs(data_hinterer["az_err"]))
data_givenpolar["az_err"] = np.minimum(np.abs(data_givenpolar["az_err"]), 360 - np.abs(data_givenpolar["az_err"]))
data_givenorient["az_err"] = np.minimum(np.abs(data_givenorient["az_err"]), 360 - np.abs(data_givenorient["az_err"]))






output_dir = './results_plots/'

# Overviews

model_names = ['Gaussian', 'Hinterer', 'Given Polar', 'Given orientation']
plot_colours = [dblue, dyellow, dgreen, dpink]
hist_colours = [dBlues, dYellows, dGreens, dPinks]

# Generate plots for each fixed inclination
for inclination in [0, 23, 45, 68, 90]:

    # ------------------------------------
    # Overview plots with background image

    # Load and display background image for overview plots
    background_image = plt.imread(f'animation_lowres/vary_inc/sim_inc{round(inclination):03d}_az0_run1.png')

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    # Loop over models
    for i, data in enumerate([data_gaussian, data_hinterer, data_givenpolar, data_givenorient]):

      data_fixed_inc = data[abs(data['inc_tru'] - inclination) <= 5]

      # Count ring vs centre
      centre_definition = data_fixed_inc["x_err"]**2 + data_fixed_inc["y_err"]**2 > 60**2
      frac_outliers = centre_definition.sum()/len(data_fixed_inc)

      mean_para = np.mean(data_fixed_inc["para_err"])
      std_para = np.std(data_fixed_inc["para_err"])
      mean_perp = np.mean(data_fixed_inc["perp_err"])
      std_perp = np.std(data_fixed_inc["perp_err"])

      ax = axs.flat[i]

      ax.scatter(data_fixed_inc["para_err"], data_fixed_inc["perp_err"], s=2, color=plot_colours[i])
      ax.axvline(mean_para, color=dred, linewidth=1, label=f'μ = ({mean_para:.2f}, {mean_perp:.2f})')
      ax.axhline(mean_perp, color=dred, linewidth=1)
      ax.add_patch(
        plt.Rectangle(
            (mean_para - std_para, mean_perp - std_perp),  # Bottom-left corner
            2 * std_para, # Width
            2 * std_perp, # Height
            color=dred,
            fill=False,
            linestyle='dashed',
            linewidth=1,
            label=f'σ = ({std_para:.2f}, {std_perp:.2f})'
        )
      )

      ax.set_facecolor('black')
      ax.imshow(background_image, cmap='gray', extent=(-500, 500, -500, 500), aspect='auto', zorder=-1)
      ax.set_xlim(-700, 700)
      ax.set_ylim(-700, 700)
      ax.set_xlabel('Δ$\parallel$')
      ax.set_ylabel('Δ$\perp$')
      ax.set_title(f"{model_names[i]}\nOverview of localisation residuals\n{np.round(frac_outliers*100, 2)}% in ring")
      ax.legend(loc='upper right')
   

    plt.tight_layout()
        
#    plt.show()
    output_filename = f"overview_{model_names[i]}_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}")
    plt.close()     

    # ------------------------------------
    # 2D scatter - centre only

    # Load and display background image for overview plots

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    # Loop over models
    for i, data in enumerate([data_gaussian, data_hinterer, data_givenpolar, data_givenorient]):

      data_fixed_inc = data[abs(data['inc_tru'] - inclination) <= 5]

      # Count ring vs centre
      centre_definition = (data_fixed_inc["x_err"]**2 + data_fixed_inc["y_err"]**2) < 60**2
      frac_outliers = centre_definition.sum()/len(data_fixed_inc)

      # Take centre blobs only
      data_centre = data_fixed_inc[centre_definition]

      mean_para = np.mean(data_centre["para_err"])
      std_para = np.std(data_centre["para_err"])
      mean_perp = np.mean(data_centre["perp_err"])
      std_perp = np.std(data_centre["perp_err"])

      ax = axs.flat[i]

      ax.scatter(data_centre["para_err"], data_centre["perp_err"], color=plot_colours[i], alpha=0.05, s=2)
      ax.axvline(mean_para, color=dred, linewidth=1, label=f'μ = ({mean_para:.2f}, {mean_perp:.2f})')
      ax.axhline(mean_perp, color=dred, linewidth=1)
      ax.add_patch(
          plt.Rectangle(
              (mean_para - std_para, mean_perp - std_perp),  # Bottom-left corner
              2 * std_para, # Width
              2 * std_perp, # Height
              color=dred,
              fill=False,
              linestyle='dashed',
              linewidth=1,
              label=f'σ = ({std_para:.2f}, {std_perp:.2f})'
          )
      )
      ax.set_xlim(-50, 50)
      ax.set_ylim(-50, 50)
      ax.set_xlabel('Δ$\parallel$')
      ax.set_ylabel('Δ$\perp$')
      ax.set_title('Localisation residuals histogram - centre only')
      ax.set_title(f"{model_names[i]}\nLocalisation residuals histogram\nCentre only")
      ax.legend(loc='upper right')

    plt.tight_layout()
        
#    plt.show()
    output_filename = f"overview_centre_{model_names[i]}_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}")
    plt.close()     

    # ------------------------------------
    # 2D scatter - para err

    # Load and display background image for overview plots

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    # Loop over models
    for i, data in enumerate([data_gaussian, data_hinterer, data_givenpolar, data_givenorient]):

      data_fixed_inc = data[abs(data['inc_tru'] - inclination) <= 5]

      # Count ring vs centre
      centre_definition = (data_fixed_inc["x_err"]**2 + data_fixed_inc["y_err"]**2) < 60**2
      frac_outliers = centre_definition.sum()/len(data_fixed_inc)

      # Take centre blobs only
      data_centre = data_fixed_inc[centre_definition]

      mean_para = np.mean(data_centre["para_err"])
      std_para = np.std(data_centre["para_err"])

      ax = axs.flat[i]

      ax.scatter(data_centre["az_tru"], data_centre["para_err"], s=10, alpha=0.1, color=plot_colours[i])
      ax.axhline(0, color=dred, linewidth=2)
      ax.set_ylim(-50, 50)
      ax.set_xlabel('φ, °')
      ax.set_ylabel('Δ$\parallel$, nm')
      ax.set_title(f"{model_names[i]}\nLocalisation residuals, parallel\nCentre only")

    plt.tight_layout()
        
#    plt.show()
    output_filename = f"para_err_{model_names[i]}_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}")
    plt.close()     

    # ------------------------------------
    # 2D scatter - perp err

    # Load and display background image for overview plots

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    # Loop over models
    for i, data in enumerate([data_gaussian, data_hinterer, data_givenpolar, data_givenorient]):

      data_fixed_inc = data[abs(data['inc_tru'] - inclination) <= 5]

      # Count ring vs centre
      centre_definition = (data_fixed_inc["x_err"]**2 + data_fixed_inc["y_err"]**2) < 60**2
      frac_outliers = centre_definition.sum()/len(data_fixed_inc)

      # Take centre blobs only
      data_centre = data_fixed_inc[centre_definition]

      mean_perp = np.mean(data_centre["perp_err"])
      std_perp = np.std(data_centre["perp_err"])

      ax = axs.flat[i]

      ax.scatter(data_centre["az_tru"], data_centre["perp_err"], s=10, alpha=0.1, color=plot_colours[i])
      ax.axhline(0, color=dred, linewidth=2)
      ax.set_ylim(-50, 50)
      ax.set_xlabel('φ, °')
      ax.set_ylabel('Δ$\perp$, nm')
      ax.set_title(f"{model_names[i]}\nLocalisation residuals, perpendicular\nCentre only")

    plt.tight_layout()
        
#    plt.show()
    output_filename = f"perp_err_{model_names[i]}_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}")
    plt.close()     

    # ------------------------------------
    # 2D scatter - az err

    # Load and display background image for overview plots

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    # Loop over models
    for i, data in enumerate([data_gaussian, data_hinterer, data_givenpolar, data_givenorient]):

      data_fixed_inc = data[abs(data['inc_tru'] - inclination) <= 5]

      # Count ring vs centre
      centre_definition = (data_fixed_inc["x_err"]**2 + data_fixed_inc["y_err"]**2) < 60**2
      frac_outliers = centre_definition.sum()/len(data_fixed_inc)

      # Take centre blobs only
      data_centre = data_fixed_inc[centre_definition]

      mean_az = np.mean(data_centre["az_err"])
      std_az = np.std(data_centre["az_err"])

      ax = axs.flat[i]

      ax.scatter(data_centre["az_tru"], data_centre["az_err"], s=10, alpha=0.1, color=plot_colours[i])
      ax.axhline(0, color=dred, linewidth=2)
      ax.set_xlabel('φ, °')
      ax.set_ylabel('Δ$\phi$, nm')
      ax.set_title(f"{model_names[i]}\nLocalisation residuals, azimuth\nCentre only")

    plt.tight_layout()
        
#    plt.show()
    output_filename = f"az_err_{model_names[i]}_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}")
    plt.close()     


'''

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
'''
