import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns


def compute_para_perp(x_errors, y_errors, azimuths):
    """Compute parallel and perpendicular errors"""
    para_errors = x_errors*np.cos(azimuths) - y_errors*np.sin(azimuths)
    perp_errors = x_errors*np.sin(azimuths) + y_errors*np.cos(azimuths)
    return para_errors, perp_errors


def compute_para_perp_std(x_stds, y_stds, azimuths):
    """Compute standard deviations for parallel and perpendicular components
    using error propagation formula for uncorrelated variables"""
    cos_az = np.cos(azimuths)
    sin_az = np.sin(azimuths)
    
    # Variance of para = (cos(az))^2 * var_x + (sin(az))^2 * var_y
    para_var = (cos_az**2) * (x_stds**2) + (sin_az**2) * (y_stds**2)
    
    # Variance of perp = (sin(az))^2 * var_x + (cos(az))^2 * var_y
    perp_var = (sin_az**2) * (x_stds**2) + (cos_az**2) * (y_stds**2)
    
    return np.sqrt(para_var), np.sqrt(perp_var)


# Get default matplotlib colours
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
dblue = default_colors[0]
dyellow = default_colors[1]
dgreen = default_colors[2]
dred = default_colors[3]
dBlues = LinearSegmentedColormap.from_list('dblue_to_white', [(1, 1, 1), dblue], N=100)
dYellows = LinearSegmentedColormap.from_list('dyellow_to_white', [(1, 1, 1), dyellow], N=100)

# Define file paths for CSV files
file_paths = [
    'fitting_results_hinterer_test_first100.csv',
    'fitting_results_hinterer_test_first100.csv',
    'fitting_results_hinterer.csv',
]

model_names = ['blah1', 'blah2', 'Hinterer']
datasets = []

# Load and process data from CSV files
subdirectory = "./"

for i, file_path in enumerate(file_paths):
    full_path = os.path.join(subdirectory, file_path)
    
    if os.path.exists(full_path):
        # Load data directly from CSV
        df = pd.read_csv(full_path)
        
        # Calculate errors if they're not already in the CSV
        if 'para_err' not in df.columns or 'perp_err' not in df.columns:
            # Convert azimuths to radians for calculation if they're in degrees
            az_rad = df['az_tru'] * np.pi / 180 if df['az_tru'].max() > 6.28 else df['az_tru']
            
            # Calculate parallel and perpendicular errors
            para_errors, perp_errors = compute_para_perp(
                df['x_tru'] - df['x_est'], 
                df['y_tru'] - df['y_est'], 
                az_rad
            )
            
            df['para_err'] = para_errors
            df['perp_err'] = perp_errors
            
            # Calculate standard deviations for para and perp if standard deviations are available
            if 'x_std' in df.columns and 'y_std' in df.columns:
                para_stds, perp_stds = compute_para_perp_std(
                    df['x_std'], 
                    df['y_std'], 
                    az_rad
                )
                
                df['para_std'] = para_stds
                df['perp_std'] = perp_stds
        
        datasets.append(df)
    else:
        print(f"Warning: File not found: {full_path}")
        # Add an empty DataFrame to maintain index alignment
        datasets.append(pd.DataFrame())

for dataset in datasets:
    # Convert rad to deg
    angle_columns = ["inc_tru", "az_tru", "inc_est", "az_est", "inc_err", "az_err"]
    dataset[angle_columns] = dataset[angle_columns] * 180 / np.pi
    
    # Also convert standard deviations if present
    if 'inc_std' in dataset.columns and 'az_std' in dataset.columns:
        std_angle_columns = ["inc_std", "az_std"]
        dataset[std_angle_columns] = dataset[std_angle_columns] * 180 / np.pi

    # Account for wrapping of inc around 180 degrees:
    dataset["inc_err"] = np.mod(dataset["inc_err"], 180)
    dataset["inc_err"] = np.minimum(np.abs(dataset["inc_err"]), 180 - np.abs(dataset["inc_err"]))

    # Account for wrapping of az around 360 degrees:
    dataset["az_err"] = np.mod(dataset["az_err"], 360)
    dataset["az_err"] = np.minimum(np.abs(dataset["az_err"]), 360 - np.abs(dataset["az_err"]))


output_dir = './'

# Generate plots for each fixed inclination
for inclination in [0, 23, 45, 68, 90]:
   
    fig, axs = plt.subplots(3, 6, figsize=(40, 15))

    for i, dataset in enumerate(datasets):

        dataset_inc = dataset[abs(dataset['inc_tru'] - inclination) <= 5]

        if dataset_inc.empty:
          continue

        IQR_multiplier = 100

        # --- First column: Para difference vs azimuth ---
        data = dataset_inc["para_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_az = dataset_inc["az_tru"][(data >= lower_bound) & (data <= upper_bound)]
        
        # Get standard deviations if available
        if 'para_std' in dataset_inc.columns:
            filtered_std = dataset_inc["para_std"][(data >= lower_bound) & (data <= upper_bound)]
            has_std = True
        else:
            has_std = False
        
        mean = np.mean(filtered_data)
        std = np.std(filtered_data, ddof=1)
        se = std / np.sqrt(np.size(filtered_data))
        
        # Main scatter plot
        axs[i, 0].scatter(filtered_az, filtered_data, s=2, alpha=0.5)
        
        # Add error bars if we have standard deviation data
        if has_std:
            # To avoid too many error bars, sample a subset of points
            n_points = len(filtered_az)
            if n_points > 100000000000:
                sample_indices = np.random.choice(n_points, 100, replace=False)
                sampled_az = filtered_az.iloc[sample_indices]
                sampled_data = filtered_data.iloc[sample_indices]
                sampled_std = filtered_std.iloc[sample_indices]
                
                # Add error bars for sampled points
                axs[i, 0].errorbar(sampled_az, sampled_data, yerr=sampled_std, 
                                 fmt='none', ecolor=dblue, alpha=0.3, capsize=2)
            else:
                # Add error bars for all points
                axs[i, 0].errorbar(filtered_az, filtered_data, yerr=filtered_std, 
                                 fmt='none', ecolor=dblue, alpha=0.3, capsize=2)
                
        axs[i, 0].axhline(mean, color=dred, linewidth=1)
        axs[i, 0].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 0].set_ylim(-20, 20)
        axs[i, 0].set_xlabel('$\phi$, °')
        axs[i, 0].set_ylabel('Δ$\parallel$, nm')
        axs[i, 0].set_title(
            f"{model_names[i]}\nParallel localisation residuals\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )

        # --- Second column: Perp difference vs azimuth ---
        data = dataset_inc["perp_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_az = dataset_inc["az_tru"][(data >= lower_bound) & (data <= upper_bound)]
        
        # Get standard deviations if available
        if 'perp_std' in dataset_inc.columns:
            filtered_std = dataset_inc["perp_std"][(data >= lower_bound) & (data <= upper_bound)]
            has_std = True
        else:
            has_std = False
            
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        
        # Main scatter plot
        axs[i, 1].scatter(filtered_az, filtered_data, s=2, alpha=0.5)
        
        # Add error bars if we have standard deviation data
        if has_std:
            # To avoid too many error bars, sample a subset of points
            n_points = len(filtered_az)
            if n_points > 100000000000:
                sample_indices = np.random.choice(n_points, 100, replace=False)
                sampled_az = filtered_az.iloc[sample_indices]
                sampled_data = filtered_data.iloc[sample_indices]
                sampled_std = filtered_std.iloc[sample_indices]
                
                # Add error bars for sampled points
                axs[i, 1].errorbar(sampled_az, sampled_data, yerr=sampled_std, 
                                 fmt='none', ecolor=dblue, alpha=0.3, capsize=2)
            else:
                # Add error bars for all points
                axs[i, 1].errorbar(filtered_az, filtered_data, yerr=filtered_std, 
                                 fmt='none', ecolor=dblue, alpha=0.3, capsize=2)
                
        axs[i, 1].axhline(mean, color=dred, linewidth=1)
        axs[i, 1].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 1].set_ylim(-20, 20)
        axs[i, 1].set_xlabel('$\phi$, °')
        axs[i, 1].set_ylabel('Δ$\perp$, nm')
        axs[i, 1].set_title(
            f"{model_names[i]}\nPerpendicular localisation residuals\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )

        # --- Third column: Inc difference vs azimuth ---
        data = dataset_inc["inc_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_az = dataset_inc["az_tru"][(data >= lower_bound) & (data <= upper_bound)]
        
        # Get standard deviations if available
        if 'inc_std' in dataset_inc.columns:
            filtered_std = dataset_inc["inc_std"][(data >= lower_bound) & (data <= upper_bound)]
            has_std = True
        else:
            has_std = False
            
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        
        # Main scatter plot
        axs[i, 2].scatter(filtered_az, filtered_data, s=2, alpha=0.5)
        
        # Add error bars if we have standard deviation data
        if has_std:
            # To avoid too many error bars, sample a subset of points
            n_points = len(filtered_az)
            if n_points > 100000000000:
                sample_indices = np.random.choice(n_points, 100, replace=False)
                sampled_az = filtered_az.iloc[sample_indices]
                sampled_data = filtered_data.iloc[sample_indices]
                sampled_std = filtered_std.iloc[sample_indices]
                
                # Add error bars for sampled points
                axs[i, 2].errorbar(sampled_az, sampled_data, yerr=sampled_std, 
                                 fmt='none', ecolor=dblue, alpha=0.3, capsize=2)
            else:
                # Add error bars for all points
                axs[i, 2].errorbar(filtered_az, filtered_data, yerr=filtered_std, 
                                 fmt='none', ecolor=dblue, alpha=0.3, capsize=2)
                
        axs[i, 2].axhline(mean, color=dred, linewidth=1)
        axs[i, 2].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 2].set_ylim(0, 90)
        axs[i, 2].set_xlabel('$\phi$, °')
        axs[i, 2].set_ylabel('Δθ, °')
        axs[i, 2].set_title(
            f"{model_names[i]}\nPolar residuals\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )

        # --- Fourth column: Az difference vs azimuth ---
        data = dataset_inc["az_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_az = dataset_inc["az_tru"][(data >= lower_bound) & (data <= upper_bound)]
        
        # Get standard deviations if available
        if 'az_std' in dataset_inc.columns:
            filtered_std = dataset_inc["az_std"][(data >= lower_bound) & (data <= upper_bound)]
            has_std = True
        else:
            has_std = False
            
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        
        # Main scatter plot
        axs[i, 3].scatter(filtered_az, filtered_data, s=2, alpha=0.5)
        
        # Add error bars if we have standard deviation data
        if has_std:
            # To avoid too many error bars, sample a subset of points
            n_points = len(filtered_az)
            if n_points > 100000000000:
                sample_indices = np.random.choice(n_points, 100, replace=False)
                sampled_az = filtered_az.iloc[sample_indices]
                sampled_data = filtered_data.iloc[sample_indices]
                sampled_std = filtered_std.iloc[sample_indices]
                
                # Add error bars for sampled points
                axs[i, 3].errorbar(sampled_az, sampled_data, yerr=sampled_std, 
                                 fmt='none', ecolor=dblue, alpha=0.3, capsize=2)
            else:
                # Add error bars for all points
                axs[i, 3].errorbar(filtered_az, filtered_data, yerr=filtered_std, 
                                 fmt='none', ecolor=dblue, alpha=0.3, capsize=2)
                
        axs[i, 3].axhline(mean, color=dred, linewidth=1)
        axs[i, 3].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 3].set_ylim(0, 180)
        axs[i, 3].set_xlabel('$\phi$, °')
        axs[i, 3].set_ylabel('Δ$\phi$, °')
        axs[i, 3].set_title(
            f"{model_names[i]}\nAzimuth residuals\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )
        
        # --- Fifth column: Photon difference vs azimuth ---
        data = dataset_inc["photon_err"]
        # IQR to remove outliers
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_multiplier * IQR
        upper_bound = Q3 + IQR_multiplier * IQR
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_az = dataset_inc["az_tru"][(data >= lower_bound) & (data <= upper_bound)]
        
        # Get standard deviations if available
        if 'photon_std' in dataset_inc.columns:
            filtered_std = dataset_inc["photon_std"][(data >= lower_bound) & (data <= upper_bound)]
            has_std = True
        else:
            has_std = False
            
        mean = np.mean(filtered_data)
        std = np.std(filtered_data)
        se = std / np.sqrt(np.size(filtered_data))
        
        # Main scatter plot
        axs[i, 4].scatter(filtered_az, filtered_data, s=2, alpha=0.5)
        
        # Add error bars if we have standard deviation data
        if has_std:
            # To avoid too many error bars, sample a subset of points
            n_points = len(filtered_az)
            if n_points > 100000000000:
                sample_indices = np.random.choice(n_points, 100, replace=False)
                sampled_az = filtered_az.iloc[sample_indices]
                sampled_data = filtered_data.iloc[sample_indices]
                sampled_std = filtered_std.iloc[sample_indices]
                
                # Add error bars for sampled points
                axs[i, 4].errorbar(sampled_az, sampled_data, yerr=sampled_std, 
                                 fmt='none', ecolor=dblue, alpha=0.3, capsize=2)
            else:
                # Add error bars for all points
                axs[i, 4].errorbar(filtered_az, filtered_data, yerr=filtered_std, 
                                 fmt='none', ecolor=dblue, alpha=0.3, capsize=2)
                
        axs[i, 4].axhline(mean, color=dred, linewidth=1)
        axs[i, 4].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 4].set_ylim(-200, 1200)
        axs[i, 4].set_xlabel('$\phi$, °')
        axs[i, 4].set_ylabel('ΔN')
        axs[i, 4].set_title(
            f"{model_names[i]}\nPhoton count residuals\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )
        
        # --- Sixth column: Log-likelihood vs azimuth ---
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
        
        # Main scatter plot
        axs[i, 5].scatter(filtered_az, filtered_data, s=2, alpha=0.5)
                
        axs[i, 5].axhline(mean, color=dred, linewidth=1)
        axs[i, 5].fill_between(np.linspace(0, 360), mean-abs(std), mean+abs(std), color='red', alpha=0.1)
        axs[i, 5].set_ylim(0, 1200)
        axs[i, 5].set_xlabel('$\phi$, °')
        axs[i, 5].set_ylabel('log-likelihood')
        axs[i, 5].set_title(
            f"{model_names[i]}\nLog-likelihood\n"
            f"μ = {mean:.4f}, σ = {std:.4f}, μ/SE = {mean/se:.4f}"
        )

    plt.tight_layout()
    output_filename = f"scatter_hinterer_inc{round(inclination)}_with_uncertainty.pdf"
    plt.savefig(f"{output_dir}{output_filename}", dpi=300)
    plt.close()
