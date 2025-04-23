import os
import re
import numpy as np
from tifffile import imread, imwrite

# Function to extract theta and phi values from filename
def extract_angles(filename):
    theta_match = re.search(r'theta(\d+\.\d+)', filename)
    phi_match = re.search(r'phi(\d+\.\d+)', filename)
    
    theta = float(theta_match.group(1)) if theta_match else 0.0
    phi = float(phi_match.group(1)) if phi_match else 0.0
    
    return (theta, phi)

# Read the filenames from the text file
with open('filenames.txt', 'r') as f:
    filenames = [line.strip() for line in f.readlines()]

# Sort filenames by theta (primary) and phi (secondary)
filenames.sort(key=extract_angles)

print(f"Found {len(filenames)} files to process.")

# Create a list to collect all images in memory
# This is fine as scientific TIF images are usually small
images = []
failed_files = []

# Load all images
for i, filename in enumerate(filenames):
    try:
        img = imread(filename)
        images.append(img)
        if i % 100 == 0:
            print(f"Loaded {i}/{len(filenames)}: {filename}")
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        failed_files.append(filename)

print(f"Successfully loaded {len(images)} out of {len(filenames)} files.")
if failed_files:
    print(f"Failed to load {len(failed_files)} files.")

# If we have images, save as a multi-page TIFF
if images:
    # Convert list of images to a 3D numpy array
    stack = np.array(images)
    
    # Path for the output file
    output_path = "stack_all_simulations.tif"
    
    # Save as multi-page TIFF
    try:
        imwrite(output_path, stack, compression='zlib', photometric='minisblack')
        print(f"Successfully created TIFF stack with {len(images)} frames at {output_path}")
    except Exception as e:
        print(f"Error saving stack: {e}")
else:
    print("No images were successfully loaded. Check file paths and image format.")
