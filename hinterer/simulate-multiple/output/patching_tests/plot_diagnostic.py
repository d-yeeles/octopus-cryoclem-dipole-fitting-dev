import numpy as np
import os
import sys
import importlib.util
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import glob

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
data_1 = pd.DataFrame(columns=[
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

datasets = [data_1]
model_names = ['hinterer 4spot']
file_paths = [
    'results_4spot.py',
]

print("Starting data import...")
SPOTS_PER_IMAGE = 4

# Importing all the data from each results file
subdirectory = "./"

for i, dataset in enumerate(datasets):
    file_path = os.path.join(subdirectory, file_paths[i])
    print(f"Processing file: {file_path}")

    if os.path.exists(file_path):
        try:
            spec = importlib.util.spec_from_file_location(model_names[i], file_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            # Get parameters if they exist in the module
            frame_parameters = {}
            if hasattr(module, 'image_size_nm'):
                frame_parameters['image_size_nm'] = module.image_size_nm
            else:
                frame_parameters['image_size_nm'] = np.sqrt(SPOTS_PER_IMAGE)*2000  # Default value
                print(f"Warning: image_size_nm not found in module, using default: {frame_parameters['image_size_nm']}")

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
                    "para_err": compute_para_perp(np.array(module.x_tru) - np.array(module.x_est), 
                                                 np.array(module.y_tru) - np.array(module.y_est), 
                                                 module.az_tru)[0],
                    "perp_err": compute_para_perp(np.array(module.x_tru) - np.array(module.x_est), 
                                                 np.array(module.y_tru) - np.array(module.y_est), 
                                                 module.az_tru)[1],
                    "photon_tru": module.photon_tru,
                    "photon_est": module.photon_est,
                    "photon_err": np.array(module.photon_tru) - np.array(module.photon_est)
                })

            print(f"Data loaded successfully: {len(newdata)} rows")

            if dataset.empty:
                datasets[i] = newdata
            else:
                datasets[i] = pd.concat([dataset, newdata], ignore_index=True)
                
        except Exception as e:
            print(f"Error processing file {file_path}: {e}")
    else:
        print(f"File not found: {file_path}")

if datasets[0].empty:
    print("No data was loaded. Check file paths and data formats.")
    sys.exit(1)

print("Processing data...")

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

# Define directory and find all image files
image_dir = './sims_4spot/png_converted/'
print(f"Looking for images in: {image_dir}")
image_paths = sorted(glob.glob(os.path.join(image_dir, '*.png')))

if not image_paths:
    print(f"No PNG files found in {image_dir}")
    if os.path.exists(image_dir):
        print(f"Contents of {image_dir}:")
        for item in os.listdir(image_dir):
            print(f"  - {item}")
    else:
        print(f"Directory {image_dir} does not exist")
    sys.exit(1)

print(f"Found {len(image_paths)} images")

# Check if the number of data points is consistent with the number of images
total_spots = len(datasets[0])
expected_images = total_spots // SPOTS_PER_IMAGE
if total_spots % SPOTS_PER_IMAGE != 0:
    print(f"Warning: The number of data points ({total_spots}) is not a multiple of spots per image ({SPOTS_PER_IMAGE}).")
    print(f"Some data may not be displayed correctly.")

if expected_images != len(image_paths):
    print(f"Warning: Expected {expected_images} images based on data ({total_spots} data points / {SPOTS_PER_IMAGE} spots per image),")
    print(f"but found {len(image_paths)} images. Some data or images may not be properly matched.")

# Generate diagnostic plots
output_dir = "./diagnostic_plots"
os.makedirs(output_dir, exist_ok=True)

# Number of images to process (limited by both data and available images)
num_images_to_process = min(expected_images, len(image_paths))

for i in range(num_images_to_process):
    try:
        image_path = image_paths[i]
        print(f"Processing image {i+1}/{num_images_to_process}: {image_path}")
        
        # Calculate the start and end indices for the data slice
        start_idx = i * SPOTS_PER_IMAGE
        end_idx = start_idx + SPOTS_PER_IMAGE
        
        # Slice the dataset to get only the relevant spots for this image
        image_data = datasets[0].iloc[start_idx:end_idx]
        
        print(f"Using data points {start_idx+1} to {end_idx} for this image")
        
        # Load the image
        image = plt.imread(image_path)
        
        # Create diagnostic plot
        plt.figure(figsize=(16, 16))
        
        # Display the image with proper extent
        plt.imshow(image, cmap='gray', 
                   extent=[-frame_parameters['image_size_nm']/2, 
                           frame_parameters['image_size_nm']/2, 
                           -frame_parameters['image_size_nm']/2, 
                           frame_parameters['image_size_nm']/2])
        
        # Add data points for this specific image
        plt.scatter(image_data["x_tru"], image_data["y_tru"], c='yellow', label='True positions', alpha=0.5)
        plt.scatter(image_data["x_est"], image_data["y_est"], c='red', label='Estimated positions', alpha=0.5)
        
        # Save the figure
        output_path = os.path.join(output_dir, f"diagnostic_frame_{i+1}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Saved diagnostic plot to: {output_path}")
        
    except Exception as e:
        print(f"Error processing image {i+1}: {e}")

print("Diagnostic plot generation completed!")
