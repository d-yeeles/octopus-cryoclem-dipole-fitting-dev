import os
import random
import re
import numpy as np
import glob
import shutil
import tifffile  # For better TIFF handling

def process_images_in_directory(input_directory, output_directory, skip_missing_params=False):
    """
    Process all TIF images in the input directory by creating a larger canvas
    and randomly positioning the original image within it.
    Saves new images and updated parameter files to the output directory.
    
    Args:
        input_directory: Directory containing original TIF files
        output_directory: Directory where processed files will be saved
        skip_missing_params: If True, process images even if parameter files are missing
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    
    # Get all TIF files in the directory
    tif_files = glob.glob(os.path.join(input_directory, "*.tif"))
    print(f"Found {len(tif_files)} TIF files.")
    
    processed_count = 0
    skipped_count = 0
    no_params_count = 0
    
    for tif_file in tif_files:
        # Extract the base name without extension
        full_base_name = os.path.splitext(os.path.basename(tif_file))[0]
        
        # Extract the theta_phi_run part from filenames like "sim_theta068_phi032_run1.tif"
        # to match with parameter files like "params_theta068_phi032_run1.m"
        if full_base_name.startswith("sim_"):
            param_base_name = full_base_name[4:]  # Remove "sim_" prefix
        else:
            param_base_name = full_base_name
        
        # Construct the parameter file path
        param_file = os.path.join(input_directory, f"params_{param_base_name}.m")
        
        # Output file paths
        output_tif_path = os.path.join(output_directory, os.path.basename(tif_file))
        output_param_path = os.path.join(output_directory, f"params_{param_base_name}.m")
        
        # Check if parameter file exists
        if not os.path.exists(param_file):
            print(f"Warning: Parameter file not found for {tif_file}.")
            print(f"  - Expected parameter file: {param_file}")
            
            if skip_missing_params:
                print(f"  - Processing image without parameter file (--skip-missing-params option)")
                process_image_only(tif_file, output_tif_path)
                processed_count += 1
                no_params_count += 1
            else:
                print(f"  - Skipping {tif_file}")
                skipped_count += 1
            
            continue
        
        # Process the image and its parameter file
        process_image_and_params(tif_file, param_file, output_tif_path, output_param_path)
        processed_count += 1
    
    print(f"Processing summary: {processed_count} files processed ({no_params_count} without parameter files), {skipped_count} files skipped.")

def process_image_and_params(input_tif_path, input_param_path, output_tif_path, output_param_path):
    """
    Process a single TIF image and its parameter file, saving results to new files.
    """
    print(f"Processing: {input_tif_path}")
    
    # Open the image using tifffile
    img_array = tifffile.imread(input_tif_path)
    
    # Get original dimensions
    orig_height, orig_width = img_array.shape[:2]
    
    # Create a new blank canvas (155x155 pixels)
    new_width, new_height = 155, 155
    
    # Ensure the original image fits within the new canvas
    if orig_width > new_width or orig_height > new_height:
        print(f"Warning: Original image {input_tif_path} is larger than the target size. Skipping.")
        return
    
    # Calculate the maximum possible offset
    max_x_offset = new_width - orig_width
    max_y_offset = new_height - orig_height
    
    # Generate random offsets
    x_offset = random.randint(0, max_x_offset)
    y_offset = random.randint(0, max_y_offset)
    
    # Create new blank image (black background)
    if len(img_array.shape) == 2:  # Grayscale image
        new_img_array = np.zeros((new_height, new_width), dtype=img_array.dtype)
    else:  # Color image
        new_img_array = np.zeros((new_height, new_width, img_array.shape[2]), dtype=img_array.dtype)
    
    # Place the original image into the new canvas at the random offset
    new_img_array[y_offset:y_offset+orig_height, x_offset:x_offset+orig_width] = img_array
    
    # Save using tifffile
    tifffile.imwrite(
        output_tif_path,
        new_img_array,
        compression='lzw',  # LZW compression as in the MATLAB code
        software='MATLAB',  # Software tag
        photometric='minisblack'  # MinIsBlack (for grayscale)
    )
    
    # Update the parameter file with the new X and Y coordinates
    update_parameter_file(input_param_path, output_param_path, x_offset, y_offset)
    
    print(f"  Created image with offsets: x={x_offset}, y={y_offset}")
    print(f"  Saved to: {output_tif_path}")

def process_image_only(input_tif_path, output_tif_path):
    """
    Process a single TIF image without an associated parameter file.
    """
    print(f"Processing image only: {input_tif_path}")
    
    # Open the image using tifffile
    img_array = tifffile.imread(input_tif_path)
    
    # Get original dimensions
    orig_height, orig_width = img_array.shape[:2]
    
    # Create a new blank canvas (155x155 pixels)
    new_width, new_height = 155, 155
    
    # Ensure the original image fits within the new canvas
    if orig_width > new_width or orig_height > new_height:
        print(f"Warning: Original image {input_tif_path} is larger than the target size. Skipping.")
        return
    
    # Calculate the maximum possible offset
    max_x_offset = new_width - orig_width
    max_y_offset = new_height - orig_height
    
    # Generate random offsets
    x_offset = random.randint(0, max_x_offset)
    y_offset = random.randint(0, max_y_offset)
    
    # Create new blank image (black background)
    if len(img_array.shape) == 2:  # Grayscale image
        new_img_array = np.zeros((new_height, new_width), dtype=img_array.dtype)
    else:  # Color image
        new_img_array = np.zeros((new_height, new_width, img_array.shape[2]), dtype=img_array.dtype)
    
    # Place the original image into the new canvas at the random offset
    new_img_array[y_offset:y_offset+orig_height, x_offset:x_offset+orig_width] = img_array
    
    # Save using tifffile
    tifffile.imwrite(
        output_tif_path,
        new_img_array,
        compression='lzw',  # LZW compression as in the MATLAB code
        software='MATLAB',  # Software tag
        photometric='minisblack'  # MinIsBlack (for grayscale)
    )
    
    print(f"  Created image with offsets: x={x_offset}, y={y_offset}")
    print(f"  Saved to: {output_tif_path} (no parameter file updated)")

def update_parameter_file(input_param_path, output_param_path, x_offset, y_offset):
    """
    Update the parameter file to account for the offset and new image size.
    
    Args:
        input_param_path: Path to the original parameter file
        output_param_path: Path to save the updated parameter file
        x_offset: The x offset of the original image in the new canvas
        y_offset: The y offset of the original image in the new canvas
    """
    # Read the parameter file
    with open(input_param_path, 'r') as file:
        content = file.read()
    
    # Parse key parameters from the file
    original_size_px_match = re.search(r'image_size_px\s*=\s*(\d+)', content)
    if not original_size_px_match:
        print(f"Warning: Could not find image_size_px in {input_param_path}")
        return False
    
    original_size_px = int(original_size_px_match.group(1))
    new_size_px = 155  # The new canvas size
    
    # Calculate the center offset (to convert from top-left origin to center origin)
    # In the original image, (0,0) is at the center of the image
    # In the new image, we need to find the new center coordinates
    original_center_x = original_size_px // 2
    original_center_y = original_size_px // 2
    
    new_center_x = new_size_px // 2
    new_center_y = new_size_px // 2
    
    # Calculate the offset from the original center to the new center
    # Taking into account the top-left positioning offset
    center_offset_x = new_center_x - (original_center_x + x_offset)
    center_offset_y = new_center_y - (original_center_y + y_offset)
    
    # Update the X and Y position arrays
    # Note: In the parameter file, positive Y is upward, so we negate the Y offset
    x_pattern = re.compile(r'(positionX_nm_array\s*=\s*\[)([^]]+)(\s*,\s*\])')
    y_pattern = re.compile(r'(positionY_nm_array\s*=\s*\[)([^]]+)(\s*,\s*\])')
    
    # Function to update position values
    def update_positions(match, offset_nm):
        prefix = match.group(1)
        values_str = match.group(2)
        suffix = match.group(3)
        
        values = [float(x.strip()) for x in values_str.split(',') if x.strip()]
        updated_values = [v - offset_nm for v in values]  # Subtract offset because origin is at center
        
        return f"{prefix}{' , '.join([str(v) for v in updated_values])}{suffix}"
    
    # Parse the pixel_size_nm to convert pixel offsets to nm
    pixel_size_match = re.search(r'pixel_size_nm\s*=\s*(\d+\.\d+)', content)
    if pixel_size_match:
        pixel_size_nm = float(pixel_size_match.group(1))
        
        # Convert pixel offset to nm offset (remember center origin)
        center_offset_x_nm = center_offset_x * pixel_size_nm
        center_offset_y_nm = center_offset_y * pixel_size_nm
        
        # Apply the updates to position arrays
        updated_content = content
        updated_content = x_pattern.sub(lambda m: update_positions(m, center_offset_x_nm), updated_content)
        updated_content = y_pattern.sub(lambda m: update_positions(m, center_offset_y_nm), updated_content)
        
        # Update image size in pixels
        updated_content = re.sub(r'(image_size_px\s*=\s*)(\d+)', 
                                lambda m: m.group(1) + str(new_size_px), 
                                updated_content)
        
        # Update image size in nm (proportionally)
        if 'image_size_nm' in updated_content:
            original_size_nm_match = re.search(r'image_size_nm\s*=\s*(\d+\.\d+)', updated_content)
            if original_size_nm_match:
                original_size_nm = float(original_size_nm_match.group(1))
                new_size_nm = original_size_nm * (new_size_px / original_size_px)
                updated_content = re.sub(r'(image_size_nm\s*=\s*)(\d+\.\d+)', 
                                       lambda m: m.group(1) + f"{new_size_nm:.3f}", 
                                       updated_content)
        
        # Write the updated content to the new file
        with open(output_param_path, 'w') as file:
            file.write(updated_content)
        
        print(f"  Updated parameter file with offsets: x_nm={center_offset_x_nm:.2f}, y_nm={center_offset_y_nm:.2f}")
        return True
    else:
        print(f"Warning: Could not find pixel_size_nm in {input_param_path}")
        return False

if __name__ == "__main__":
    import sys
    
    # Check for --skip-missing-params flag
    skip_missing_params = False
    if "--skip-missing-params" in sys.argv:
        skip_missing_params = True
        sys.argv.remove("--skip-missing-params")
    
    # Simple command line argument handling to match your usage pattern
    if len(sys.argv) < 2:
        print("Usage: python make_bigger_all_3.py input_directory [output_directory] [--skip-missing-params]")
        print("  --skip-missing-params: Process images even if their parameter files are missing")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    
    # If output directory is provided, use it; otherwise, add "_processed" to input directory
    if len(sys.argv) >= 3:
        output_dir = sys.argv[2]
    else:
        output_dir = input_dir + "_processed"
    
    print(f"Processing TIF images from: {input_dir}")
    print(f"Saving results to: {output_dir}")
    if skip_missing_params:
        print("Option: Will process images even if parameter files are missing")
    
    process_images_in_directory(input_dir, output_dir, skip_missing_params)
    print("Processing complete!")
