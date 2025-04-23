#!/usr/bin/env python3
"""
Headless Dipole Position Visualizer

This version saves directly to file without interactive display to avoid segmentation faults.

Usage:
    python position_visualizer_headless.py <path_to_tif_file> <path_to_params_file> <output_file>

Example:
    python position_visualizer_headless.py sim_theta068_phi000_run1.tif params_theta068_phi000_run1.m visualization.png

Requirements:
    pip install numpy matplotlib pillow tifffile
"""

import sys
import os
import re
import numpy as np
import matplotlib
# Use non-interactive Agg backend to avoid display issues
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Try to import tifffile
try:
    import tifffile
    HAVE_TIFFFILE = True
except ImportError:
    HAVE_TIFFFILE = False
    print("Warning: tifffile module not found. Using matplotlib.pyplot.imread instead.")

def extract_array_values(line):
    """Extract numerical values from a line containing an array definition."""
    match = re.search(r'\[(.*)\]', line)
    if match:
        values_str = match.group(1)
        values = [float(x.strip()) for x in values_str.split(',') if x.strip()]
        return values
    return []

def parse_params_file(params_file):
    """Parse the MATLAB params file to extract positions and image parameters."""
    params = {}
    
    with open(params_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('%') or not line:
                continue
            
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                
                if '[' in value and ']' in value:
                    params[key] = extract_array_values(value)
                else:
                    try:
                        params[key] = float(value.strip())
                    except ValueError:
                        params[key] = value.strip()
    
    return params

def read_image(file_path):
    """Read an image file with appropriate method based on available libraries."""
    if HAVE_TIFFFILE and file_path.lower().endswith(('.tif', '.tiff')):
        print("Reading with tifffile...")
        try:
            return tifffile.imread(file_path)
        except Exception as e:
            print(f"Error with tifffile: {e}")
    
    print("Reading with matplotlib...")
    try:
        return plt.imread(file_path)
    except Exception as e:
        print(f"Error reading image: {e}")
        raise e

def normalize_image(image):
    """Normalize image to 0-1 range for display."""
    img_float = image.astype(np.float64)
    min_val = np.min(img_float)
    max_val = np.max(img_float)
    
    if max_val == min_val:
        return np.zeros_like(img_float)
    
    return (img_float - min_val) / (max_val - min_val)

def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <tif_file> <params_file> <output_file>")
        sys.exit(1)
    
    tif_file = sys.argv[1]
    params_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Load and normalize image
    try:
        image = read_image(tif_file)
        print(f"Image loaded: shape={image.shape}, dtype={image.dtype}")
        print(f"Min={np.min(image)}, Max={np.max(image)}")
        
        normalized_image = normalize_image(image)
        
        # Check if image is mostly black
        black_percentage = np.sum(normalized_image < 0.02) / normalized_image.size * 100
        if black_percentage > 95:
            print(f"Warning: Image is {black_percentage:.1f}% black. Trying log normalization...")
            
            # Try log normalization
            eps = 1e-10
            log_image = np.log(image.astype(np.float64) + eps)
            normalized_image = normalize_image(log_image)
    except Exception as e:
        print(f"Failed to load or process image: {e}")
        print("Creating a blank image for visualization...")
        # Create blank image
        image = np.zeros((155, 155), dtype=np.float32)
        normalized_image = image.copy()
    
    # Parse params file
    try:
        params = parse_params_file(params_file)
        
        # Extract key parameters
        pixel_size_nm = float(params.get('pixel_size_nm', 52))
        image_size_px = int(params.get('image_size_px', 155))
        
        # Get position arrays
        positions_x_nm = params.get('positionX_nm_array', [])
        positions_y_nm = params.get('positionY_nm_array', [])
        
        print(f"Parsed {len(positions_x_nm)} dipole positions")
        print(f"Image parameters: size={image_size_px}px, pixel_size={pixel_size_nm}nm")
    except Exception as e:
        print(f"Failed to parse params file: {e}")
        sys.exit(1)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Display image
    ax.imshow(normalized_image, cmap='viridis', origin='lower')
    
    # Calculate center of image (in pixels)
    center_px = image_size_px // 2
    
    print(f"Image center: pixel ({center_px}, {center_px})")
    
    # Plot each dipole position
    for i, (x_nm, y_nm) in enumerate(zip(positions_x_nm, positions_y_nm)):
        # Convert from nm to pixels
        x_px = x_nm / pixel_size_nm + center_px
        y_px = center_px - y_nm / pixel_size_nm  # Flip y-axis
        
        # Plot position
        circle = Circle((x_px, y_px), radius=5, color='red', fill=True, alpha=0.7)
        ax.add_patch(circle)
        
        # Annotate with position index
        ax.text(x_px + 7, y_px, f"{i+1}", color='white', fontsize=12, 
                bbox=dict(facecolor='black', alpha=0.7))
        
        # Print coordinates for verification
        print(f"Dipole {i+1}: nm({x_nm:.2f}, {y_nm:.2f}) -> px({x_px:.2f}, {y_px:.2f})")
    
    # Add title
    ax.set_title(f"Dipole Positions\nImage: {os.path.basename(tif_file)}")
    
    # Add coordinate system indicators
    arrow_props = dict(arrowstyle='->', lw=2, color='cyan')
    ax.annotate('', xy=(center_px + 50, center_px), xytext=(center_px, center_px), 
                arrowprops=arrow_props)
    ax.annotate('', xy=(center_px, center_px + 50), xytext=(center_px, center_px), 
                arrowprops=arrow_props)
    ax.text(center_px + 55, center_px, 'x', color='cyan', fontsize=14)
    ax.text(center_px, center_px + 55, 'y', color='cyan', fontsize=14)
    
    # Add pixel grid for reference
    ax.grid(which='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
    
    plt.tight_layout()
    
    # Save file
    print(f"Saving visualization to {output_file}")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)  # Explicitly close figure to free memory
    print("Done!")

if __name__ == "__main__":
    main()
