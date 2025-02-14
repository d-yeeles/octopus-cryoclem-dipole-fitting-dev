import numpy as np
import os
import sys
import importlib.util
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.image as mpimg
import subprocess
from sklearn.cluster import KMeans
from scipy.spatial import ConvexHull
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

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

    # Count ring vs centre
    centre_definition_gaussian = data_gaussian_fixed_inc["x_err"]**2 + data_gaussian_fixed_inc["y_err"]**2 > 60**2
    frac_outliers_gaussian = centre_definition_gaussian.sum()/len(data_gaussian_fixed_inc)
    centre_definition_hinterer = data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 > 60**2
    frac_outliers_hinterer = centre_definition_hinterer.sum()/len(data_hinterer_fixed_inc)

    fig, axs = plt.subplots(2, 1, figsize=(8, 16))

    # Gaussian - Clustering and visualization
    para_err_centre = data_gaussian_fixed_inc["para_err"][data_gaussian_fixed_inc["x_err"]**2 + data_gaussian_fixed_inc["y_err"]**2 < 60**2]
    perp_err_centre = data_gaussian_fixed_inc["perp_err"][data_gaussian_fixed_inc["x_err"]**2 + data_gaussian_fixed_inc["y_err"]**2 < 60**2]

    # Perform DBSCAN clustering (adjust eps to make smaller clusters)
    data_points = np.column_stack((para_err_centre, perp_err_centre))
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data_points)
    dbscan = DBSCAN(eps=0.23, min_samples=100)  # Smaller eps and higher min_samples for smaller clusters
    labels = dbscan.fit_predict(data_scaled)

    # Plot Gaussian data with clusters and convex hulls
    unique_labels = set(labels)
    for label in unique_labels:
        cluster_points = data_points[labels == label]
        if label != -1:  # Skip noise points
            axs[0].scatter(cluster_points[:, 0], cluster_points[:, 1])

            # Draw convex hull boundary for each cluster
            if len(cluster_points) > 2:  # Convex hull requires at least 3 points
                hull = ConvexHull(cluster_points)
                for simplex in hull.simplices:
                    axs[0].plot(cluster_points[simplex, 0], cluster_points[simplex, 1], 'k-', lw=2)

            # Add the number of members in the cluster
            axs[0].text(np.mean(cluster_points[:, 0]), np.mean(cluster_points[:, 1]),
                        f'{len(cluster_points)}', fontsize=12, ha='center', color='black')

    # Plot Gaussian - centre only histogram
    c = axs[0].hist2d(para_err_centre, perp_err_centre, bins=100, cmap=dBlues)
    cbar = fig.colorbar(c[3], ax=axs[0], label='Counts')
    axs[0].set_xlim(-50, 50)
    axs[0].set_ylim(-50, 50)
    axs[0].set_xlabel('Δ$\parallel$')
    axs[0].set_ylabel('Δ$\perp$')
    axs[0].set_title('Localisation residuals histogram - Gaussian centre only')

    # Hinterer - Clustering and visualization
    para_err_centre = data_hinterer_fixed_inc["para_err"][data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 < 60**2]
    perp_err_centre = data_hinterer_fixed_inc["perp_err"][data_hinterer_fixed_inc["x_err"]**2 + data_hinterer_fixed_inc["y_err"]**2 < 60**2]

    # Perform DBSCAN clustering (adjust eps to make smaller clusters)
    data_points = np.column_stack((para_err_centre, perp_err_centre))
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data_points)
    dbscan = DBSCAN(eps=0.22, min_samples=150)  # Smaller eps and higher min_samples for smaller clusters
    labels = dbscan.fit_predict(data_scaled)

    # Plot Hinterer data with clusters and convex hulls
    unique_labels = set(labels)
    for label in unique_labels:
        cluster_points = data_points[labels == label]
        if label != -1:  # Skip noise points
            axs[1].scatter(cluster_points[:, 0], cluster_points[:, 1])

            # Draw convex hull boundary for each cluster
            if len(cluster_points) > 2:  # Convex hull requires at least 3 points
                hull = ConvexHull(cluster_points)
                for simplex in hull.simplices:
                    axs[1].plot(cluster_points[simplex, 0], cluster_points[simplex, 1], 'k-', lw=2)

            # Add the number of members in the cluster
            axs[1].text(np.mean(cluster_points[:, 0]), np.mean(cluster_points[:, 1]),
                        f'{len(cluster_points)}', fontsize=12, ha='center', color='black')

    # Plot Hinterer - centre only histogram using dyellow
    c = axs[1].hist2d(para_err_centre, perp_err_centre, bins=100, cmap=dYellows)
    cbar = fig.colorbar(c[3], ax=axs[1], label='Counts')
    axs[1].set_xlim(-50, 50)
    axs[1].set_ylim(-50, 50)
    axs[1].set_xlabel('Δ$\parallel$')
    axs[1].set_ylabel('Δ$\perp$')
    axs[1].set_title('Localisation residuals histogram - Hinterer centre only')

    plt.tight_layout()
    plt.show()

