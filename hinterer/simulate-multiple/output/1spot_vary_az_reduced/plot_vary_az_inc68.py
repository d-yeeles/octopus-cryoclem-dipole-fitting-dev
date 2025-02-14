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
    "perp_err",
    "obj_tru",
    "obj_est",
    "obj_err"
])

data_hinterer = data_gaussian.copy()
data_hinterer_ps = data_gaussian.copy()




# Importing all the data from each results file
subdirectory = "./"
for filename in os.listdir(subdirectory):
#    if filename.startswith("fitting_results_inc") and filename.endswith("_gaussian.py"):
    if filename.startswith("fitting_results_gaussian_highN_fmincon"):

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
                "perp_err": compute_para_perp(gaussian.x_err, gaussian.y_err, gaussian.az_tru)[1],
                "obj_tru": gaussian.obj_tru,
                "obj_est": gaussian.obj_est,
                "obj_err": gaussian.obj_err
            })

        if data_gaussian.empty:
            data_gaussian = newdata_gaussian
        else:
            #data_gaussian = data_gaussian.reset_index(drop=True)
            data_gaussian = pd.concat([data_gaussian, newdata_gaussian], ignore_index=True)



subdirectory = "./"
for filename in os.listdir(subdirectory):
    if filename.startswith("fitting_results_hinterer_highN_swarm_inc68only."):

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


subdirectory = "./"
for filename in os.listdir(subdirectory):
    if filename.startswith("fitting_results_hinterer_highN_swarm_reparam"):

        file_path = os.path.join(subdirectory, filename)
        
        module_name = "hinterer_ps"
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        hinterer_ps = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(hinterer_ps)
        
        newdata_hinterer_ps = pd.DataFrame({
                "x_tru": hinterer_ps.x_tru,
                "y_tru": hinterer_ps.y_tru,
                "inc_tru": hinterer_ps.inc_tru,
                "az_tru": hinterer_ps.az_tru,
                "x_est": hinterer_ps.x_est,
                "y_est": hinterer_ps.y_est,
                "inc_est": hinterer_ps.inc_est,
                "az_est": hinterer_ps.az_est,
                "x_err": hinterer_ps.x_err,
                "y_err": hinterer_ps.y_err,
                "inc_err": hinterer_ps.inc_err,
                "az_err": hinterer_ps.az_err,
                "para_err": compute_para_perp(hinterer_ps.x_err, hinterer_ps.y_err, hinterer_ps.az_tru)[0],
                "perp_err": compute_para_perp(hinterer_ps.x_err, hinterer_ps.y_err, hinterer_ps.az_tru)[1],
                "obj_tru": hinterer_ps.obj_tru,
                "obj_est": hinterer_ps.obj_est,
                "obj_err": hinterer_ps.obj_err

            })

        if data_hinterer_ps.empty:
            data_hinterer_ps = newdata_hinterer_ps
        else:
            #data_hinterer_ps = data_hinterer_ps.reset_index(drop=True)
            data_hinterer_ps = pd.concat([data_hinterer_ps, newdata_hinterer_ps], ignore_index=True)



# Convert rad to deg
angle_columns = ["inc_tru", "az_tru", "inc_est", "az_est", "inc_err", "az_err"]
data_gaussian[angle_columns] = data_gaussian[angle_columns] * 180 / np.pi
data_hinterer[angle_columns] = data_hinterer[angle_columns] * 180 / np.pi
data_hinterer_ps[angle_columns] = data_hinterer_ps[angle_columns] * 180 / np.pi

# Account for wrapping of inc around 180 degrees:
data_gaussian["inc_err"] = np.mod(data_gaussian["inc_err"], 180)
data_gaussian["inc_err"] = np.minimum(np.abs(data_gaussian["inc_err"]), 180 - np.abs(data_gaussian["inc_err"]))
data_hinterer["inc_err"] = np.mod(data_hinterer["inc_err"], 180)
data_hinterer["inc_err"] = np.minimum(np.abs(data_hinterer["inc_err"]), 180 - np.abs(data_hinterer["inc_err"]))
data_hinterer_ps["inc_err"] = np.mod(data_hinterer_ps["inc_err"], 180)
data_hinterer_ps["inc_err"] = np.minimum(np.abs(data_hinterer_ps["inc_err"]), 180 - np.abs(data_hinterer_ps["inc_err"]))

# Account for wrapping of az around 360 degrees:
data_gaussian["az_err"] = np.mod(data_gaussian["az_err"], 360)
data_gaussian["az_err"] = np.minimum(np.abs(data_gaussian["az_err"]), 360 - np.abs(data_gaussian["az_err"]))
data_hinterer["az_err"] = np.mod(data_hinterer["az_err"], 360)
data_hinterer["az_err"] = np.minimum(np.abs(data_hinterer["az_err"]), 360 - np.abs(data_hinterer["az_err"]))
data_hinterer_ps["az_err"] = np.mod(data_hinterer_ps["az_err"], 360)
data_hinterer_ps["az_err"] = np.minimum(np.abs(data_hinterer_ps["az_err"]), 360 - np.abs(data_hinterer_ps["az_err"]))

output_dir = './'
# Overviews

# Generate plots for each fixed inclination
for inclination in [68]:

    data_gaussian_fixed_inc = data_gaussian[abs(data_gaussian['inc_tru'] - inclination) <= 5]
    data_hinterer_fixed_inc = data_hinterer[abs(data_hinterer['inc_tru'] - inclination) <= 5]
    data_hinterer_ps_fixed_inc = data_hinterer_ps[abs(data_hinterer_ps['inc_tru'] - inclination) <= 5]

    # Load and display background image
    background_image = plt.imread(f'../animation_lowres_highN/vary_inc/sim_inc{round(inclination):03d}_az0_run1.png')

    fig, axs = plt.subplots(3, 3, figsize=(18, 15))
    
    # Gaussian
    # Plots - full

    mean1 = np.mean(data_gaussian_fixed_inc["para_err"])
    std1 = np.std(data_gaussian_fixed_inc["para_err"])
    mean2 = np.mean(data_gaussian_fixed_inc["perp_err"])
    std2 = np.std(data_gaussian_fixed_inc["perp_err"])
    axs[0, 0].scatter(data_gaussian_fixed_inc["para_err"], data_gaussian_fixed_inc["perp_err"], s=2, color=dgreen)
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

    axs[0, 0].set_xlim(-350, 350)
    axs[0, 0].set_ylim(-350, 350)
    axs[0, 0].set_xlabel('Δ$\parallel$')
    axs[0, 0].set_ylabel('Δ$\perp$')
    axs[0, 0].set_title("Gaussian fmincon\nOverview of localisation residuals")
    axs[0, 0].legend(loc='upper right')


    # Plots - zoom 40nm

    c = axs[0, 1].scatter(data_gaussian_fixed_inc["para_err"], data_gaussian_fixed_inc["perp_err"], s=2, c=data_gaussian_fixed_inc["az_tru"], cmap='hsv')
    cbar = fig.colorbar(c)
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
    axs[0, 1].set_xlim(-40, 40)
    axs[0, 1].set_ylim(-40, 40)
    axs[0, 1].set_xlabel('Δ$\parallel$')
    axs[0, 1].set_ylabel('Δ$\perp$')
    axs[0, 1].set_title('Gaussian fmincon\nZoomed to $\pm$ 40nm')
    axs[0, 1].legend(loc='upper right')
   
    # Plots - zoom 2σ
    c = axs[0, 2].scatter(data_gaussian_fixed_inc["para_err"], data_gaussian_fixed_inc["perp_err"], s=2, c=data_gaussian_fixed_inc["az_tru"], cmap='hsv')
    cbar = fig.colorbar(c)
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
    axs[0, 2].set_xlim(mean1-2*std1, mean1+2*std1)
    axs[0, 2].set_ylim(mean2-2*std2, mean2+2*std2)
    axs[0, 2].set_xlabel('Δ$\parallel$')
    axs[0, 2].set_ylabel('Δ$\perp$')
    axs[0, 2].set_title('Gaussian fmincon\nZoomed to 2σ')
    axs[0, 2].legend(loc='upper right')





    # hinterer
    # Plots - full

    mean1 = np.mean(data_hinterer_fixed_inc["para_err"])
    std1 = np.std(data_hinterer_fixed_inc["para_err"])
    mean2 = np.mean(data_hinterer_fixed_inc["perp_err"])
    std2 = np.std(data_hinterer_fixed_inc["perp_err"])
    axs[1, 0].scatter(data_hinterer_fixed_inc["para_err"], data_hinterer_fixed_inc["perp_err"], s=2, color=dgreen)
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

    axs[1, 0].set_xlim(-350, 350)
    axs[1, 0].set_ylim(-350, 350)
    axs[1, 0].set_xlabel('Δ$\parallel$')
    axs[1, 0].set_ylabel('Δ$\perp$')
    axs[1, 0].set_title("hinterer fmincon\nOverview of localisation residuals")
    axs[1, 0].legend(loc='upper right')


    # Plots - zoom 40nm

    c = axs[1, 1].scatter(data_hinterer_fixed_inc["para_err"], data_hinterer_fixed_inc["perp_err"], s=2, c=data_hinterer_fixed_inc["az_tru"], cmap='hsv')
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
    axs[1, 1].set_xlim(-40, 40)
    axs[1, 1].set_ylim(-40, 40)
    axs[1, 1].set_xlabel('Δ$\parallel$')
    axs[1, 1].set_ylabel('Δ$\perp$')
    axs[1, 1].set_title('hinterer fmincon\nZoomed to $\pm$ 40nm')
    axs[1, 1].legend(loc='upper right')
   
    # Plots - zoom 2σ
    c = axs[1, 2].scatter(data_hinterer_fixed_inc["para_err"], data_hinterer_fixed_inc["perp_err"], s=2, c=data_hinterer_fixed_inc["az_tru"], cmap='hsv')
    cbar = fig.colorbar(c)
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
    axs[1, 2].set_xlim(mean1-2*std1, mean1+2*std1)
    axs[1, 2].set_ylim(mean2-2*std2, mean2+2*std2)
    axs[1, 2].set_xlabel('Δ$\parallel$')
    axs[1, 2].set_ylabel('Δ$\perp$')
    axs[1, 2].set_title('hinterer fmincon\nZoomed to 2σ')
    axs[1, 2].legend(loc='upper right')



    # hinterer
    # Plots - full

    mean1 = np.mean(data_hinterer_ps_fixed_inc["para_err"])
    std1 = np.std(data_hinterer_ps_fixed_inc["para_err"])
    mean2 = np.mean(data_hinterer_ps_fixed_inc["perp_err"])
    std2 = np.std(data_hinterer_ps_fixed_inc["perp_err"])
    axs[2, 0].scatter(data_hinterer_ps_fixed_inc["para_err"], data_hinterer_ps_fixed_inc["perp_err"], s=2, color=dgreen)
    axs[2, 0].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[2, 0].axhline(mean2, color=dred, linewidth=1)
    axs[2, 0].add_patch(
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

    axs[2, 0].set_facecolor('black')
    axs[2, 0].imshow(background_image, cmap='gray', extent=(-500, 500, -500, 500), aspect='auto', zorder=-1)

    axs[2, 0].set_xlim(-350, 350)
    axs[2, 0].set_ylim(-350, 350)
    axs[2, 0].set_xlabel('Δ$\parallel$')
    axs[2, 0].set_ylabel('Δ$\perp$')
    axs[2, 0].set_title("hinterer swarm\nOverview of localisation residuals")
    axs[2, 0].legend(loc='upper right')


    # Plots - zoom 40nm

    c = axs[2, 1].scatter(data_hinterer_ps_fixed_inc["para_err"], data_hinterer_ps_fixed_inc["perp_err"], s=2, c=data_hinterer_ps_fixed_inc["az_tru"], cmap='hsv')
    cbar = fig.colorbar(c)
    axs[2, 1].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[2, 1].axhline(mean2, color=dred, linewidth=1)
    axs[2, 1].add_patch(
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
    axs[2, 1].set_xlim(-40, 40)
    axs[2, 1].set_ylim(-40, 40)
    axs[2, 1].set_xlabel('Δ$\parallel$')
    axs[2, 1].set_ylabel('Δ$\perp$')
    axs[2, 1].set_title('hinterer swarm\nZoomed to $\pm$ 40nm')
    axs[2, 1].legend(loc='upper right')
   

    # Plots - zoom 2σ
    c = axs[2, 2].scatter(data_hinterer_ps_fixed_inc["para_err"], data_hinterer_ps_fixed_inc["perp_err"], s=2, c=data_hinterer_ps_fixed_inc["az_tru"], cmap='hsv')
    cbar = fig.colorbar(c)
    axs[2, 2].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
    axs[2, 2].axhline(mean2, color=dred, linewidth=1)
    axs[2, 2].add_patch(
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

#    # Adding little text labels to show az
#    for row in data_hinterer_ps_fixed_inc.itertuples(index=False):
#        axs[2, 2].text(
#            row.para_err + 0.001, 
#            row.perp_err + 0.001,
#            f'{row.az_tru:.1f}',  # Format as float with 1 decimal place
#            fontsize=6, color='black', alpha=0.7
#        )

    axs[2, 2].set_xlim(mean1-2*std1, mean1+2*std1)
    axs[2, 2].set_ylim(mean2-2*std2, mean2+2*std2)
    axs[2, 2].set_xlabel('Δ$\parallel$')
    axs[2, 2].set_ylabel('Δ$\perp$')
    axs[2, 2].set_title('hinterer swarm\nZoomed to 2σ')
    axs[2, 2].legend(loc='upper right')







    plt.tight_layout()
    #plt.show()
    output_filename = f"test_overview_inc{round(inclination)}.png"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()  




