import numpy as np
import os
import sys
import importlib.util
import matplotlib.pyplot as plt
import pandas as pd
import tifffile

# Initialize output directory
output_dir = './'

def process_tif_stack(tif_path, dataset):
    # Load the TIF stack
    print(f"Loading TIF stack: {tif_path}")
    tif_stack = tifffile.imread(tif_path)
    print(f"TIF stack shape: {tif_stack.shape}")
    
    # Determine number of frames
    num_frames = tif_stack.shape[0]
    
    # Get original image dimensions
    height_px, width_px = tif_stack.shape[1:] if len(tif_stack.shape) > 2 else tif_stack.shape
    
    # Create output directory for individual PNGs
    output_subdir = os.path.join(output_dir, 'overlay_frames')
    os.makedirs(output_subdir, exist_ok=True)
    print(f"Created output directory: {output_subdir}")
    
    # Pixel dimensions for coordinate calculations
    pixel_size = 52  # nm per pixel
    
    # Calculate extent in nanometers
    x_min = -width_px * pixel_size / 2
    x_max = width_px * pixel_size / 2
    y_min = -height_px * pixel_size / 2
    y_max = height_px * pixel_size / 2
    
    print(f"Image extent: x={x_min} to {x_max}, y={y_min} to {y_max}")
    
    # Print table of frame numbers and point counts
    for frame_num in sorted(dataset['frame'].unique()):
        count = len(dataset[dataset['frame'] == frame_num])
        print(f"Data frame {frame_num}: {count} points")
    
    # Process each frame
    for frame_idx in range(num_frames):
        # Frame index in dataset starts at 1, not 0
        data_frame_idx = frame_idx + 1
        
        # Get points for this frame
        frame_data = dataset[dataset['frame'] == data_frame_idx]
        
        print(f"Processing TIF frame {frame_idx} (data frame {data_frame_idx}): {len(frame_data)} points")
        
        # Create a figure
        dpi = 100
        fig, ax = plt.subplots(figsize=(width_px/dpi, height_px/dpi), dpi=dpi)
        
        # Display the current frame as background
        ax.imshow(tif_stack[frame_idx], cmap='gray', extent=[x_min, x_max, y_min, y_max])
        
        if len(frame_data) > 0:
            # Print the first point coordinates for verification
            print(f"  First point coordinates:")
            print(f"    True: ({frame_data['x_tru'].iloc[0]:.1f}, {frame_data['y_tru'].iloc[0]:.1f})")
            print(f"    Est:  ({frame_data['x_est'].iloc[0]:.1f}, {frame_data['y_est'].iloc[0]:.1f})")
            
            # True positions (red circles)
            ax.scatter(frame_data["x_tru"], frame_data["y_tru"], 
                     s=1, color='red')
            
            # Estimated positions (green squares)
            ax.scatter(frame_data["x_est"], frame_data["y_est"], 
                     s=1, color='lime')
            
#            # Connecting lines (yellow)
#            for _, row in frame_data.iterrows():
#                ax.plot([row["x_tru"], row["x_est"]], 
#                      [row["y_tru"], row["y_est"]], 
#                      'yellow', linewidth=1)
        
        # Set axes limits to match the image extent
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        
        # Clean up plot
        ax.set_axis_off()
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        
        # Save as PNG
        output_filename = os.path.join(output_subdir, f'frame_{(frame_idx+1):04d}.png')
        plt.savefig(output_filename, dpi=dpi, bbox_inches='tight', pad_inches=0)
        
        # Close the plot
        plt.close(fig)
    
    print(f"Created {num_frames} PNG files in {output_subdir}")

# Main script
if __name__ == "__main__":
    print("Reading data file...")
    # Load the module containing your data
    file_path = './results_real.py'
    spec = importlib.util.spec_from_file_location('results_real', file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    # Create dataset
    data = pd.DataFrame({
        "frame": module.frame,
        "x_tru": module.x_tru,
        "y_tru": module.y_tru,
        "x_est": module.x_est,
        "y_est": module.y_est,
    })
    
    print(f"Loaded {len(data)} points from data file")
    
    # Process TIF stack
    tif_path = "stack_all_simulations.tif"  # Replace with your actual file
    process_tif_stack(tif_path, data)
