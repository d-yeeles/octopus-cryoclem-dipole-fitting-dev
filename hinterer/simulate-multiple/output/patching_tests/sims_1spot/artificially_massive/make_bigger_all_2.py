#!/usr/bin/env python3
"""
Script to process all TIFF images in a directory, placing each original image
at a random position within a new 8000×8000 pixel canvas, and updating the
associated parameter files with the new feature coordinates.

Usage:
    python update_params_random_offset_tiffs.py input_dir output_dir
    
    - input_dir: directory containing TIFF images and parameter files to process
    - output_dir: directory to save the processed TIFF images and updated parameter files
"""

import os
import sys
import random
import re
from PIL import Image
import glob
import shutil

def parse_parameter_file(param_file_path):
    """
    Parse parameter file to extract position and other information.
    
    Parameters:
    -----------
    param_file_path : str
        Path to the parameter file
        
    Returns:
    --------
    dict
        Dictionary containing the parsed parameters
    """
    params = {}
    
    try:
        with open(param_file_path, 'r') as f:
            content = f.read()
            
            # Extract key parameters using regular expressions
            # Position X values
            x_match = re.search(r'positionX_nm_array\s*=\s*\[(.*?)\]', content)
            if x_match:
                x_values = [float(x.strip()) for x in x_match.group(1).split(',') if x.strip()]
                params['positionX_nm_array'] = x_values
            
            # Position Y values
            y_match = re.search(r'positionY_nm_array\s*=\s*\[(.*?)\]', content)
            if y_match:
                y_values = [float(y.strip()) for y in y_match.group(1).split(',') if y.strip()]
                params['positionY_nm_array'] = y_values
            
            # Image size in nanometers
            size_nm_match = re.search(r'image_size_nm\s*=\s*(\d+\.?\d*)', content)
            if size_nm_match:
                params['image_size_nm'] = float(size_nm_match.group(1))
            
            # Image size in pixels
            size_px_match = re.search(r'image_size_px\s*=\s*(\d+)', content)
            if size_px_match:
                params['image_size_px'] = int(size_px_match.group(1))
            
            # Pixel size in nanometers
            pixel_size_match = re.search(r'pixel_size_nm\s*=\s*(\d+\.?\d*)', content)
            if pixel_size_match:
                params['pixel_size_nm'] = float(pixel_size_match.group(1))
            
        return params
    except Exception as e:
        print(f"Error parsing parameter file {param_file_path}: {e}", file=sys.stderr)
        return None

def update_parameter_file(param_file_path, output_param_path, offset_x_nm, offset_y_nm):
    """
    Update parameter file with new position values.
    
    Parameters:
    -----------
    param_file_path : str
        Path to the original parameter file
    output_param_path : str
        Path to save the updated parameter file
    offset_x_nm : float
        X offset in nanometers to add to all X positions
    offset_y_nm : float
        Y offset in nanometers to add to all Y positions
        
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        # Read the entire original file content
        with open(param_file_path, 'r') as f:
            lines = f.readlines()
        
        # Process the file line by line
        updated_lines = []
        for line in lines:
            # Handle X positions
            if 'positionX_nm_array' in line:
                # Parse the current values
                matches = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', line)
                if matches:
                    # Update each value with the offset
                    new_values = [float(x) + offset_x_nm for x in matches]
                    # Reconstruct the line with the same format
                    new_values_str = ', '.join([f"{x:.6e}" for x in new_values]) + ','
                    new_line = f"positionX_nm_array = [{new_values_str}]\n"
                    updated_lines.append(new_line)
                else:
                    updated_lines.append(line)
            # Handle Y positions
            elif 'positionY_nm_array' in line:
                # Parse the current values
                matches = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', line)
                if matches:
                    # Update each value with the offset
                    new_values = [float(y) + offset_y_nm for y in matches]
                    # Reconstruct the line with the same format
                    new_values_str = ', '.join([f"{y:.6e}" for y in new_values]) + ','
                    new_line = f"positionY_nm_array = [{new_values_str}]\n"
                    updated_lines.append(new_line)
                else:
                    updated_lines.append(line)
            else:
                updated_lines.append(line)
        
        # Write the updated content to the output file in the output directory
        with open(output_param_path, 'w') as f:
            f.writelines(updated_lines)
        
        return True
    except Exception as e:
        print(f"Error updating parameter file {param_file_path}: {e}", file=sys.stderr)
        return False

def find_param_file_for_tiff(tiff_path):
    """
    Find the corresponding parameter file for a TIFF image.
    
    Parameters:
    -----------
    tiff_path : str
        Path to the TIFF image
        
    Returns:
    --------
    str or None
        Path to the parameter file if found, None otherwise
    """
    tiff_basename = os.path.basename(tiff_path)
    tiff_name_no_ext = os.path.splitext(tiff_basename)[0]
    
    # Handle the specific naming pattern: sim_theta068_phi308_run1.tif -> params_theta068_phi308_run1.m
    if tiff_name_no_ext.startswith("sim_"):
        param_name = f"params_{tiff_name_no_ext[4:]}.m"  # Remove 'sim_' prefix
        param_path = os.path.join(os.path.dirname(tiff_path), param_name)
        if os.path.exists(param_path):
            return param_path
    
    # Common parameter file patterns as fallbacks
    param_patterns = [
        f"params_{tiff_name_no_ext}.m",
        f"{tiff_name_no_ext}_params.m",
        f"params_{tiff_name_no_ext.replace('run', '')}.m"
    ]
    
    # Search for parameter files in the same directory
    tiff_dir = os.path.dirname(tiff_path)
    for pattern in param_patterns:
        param_path = os.path.join(tiff_dir, pattern)
        if os.path.exists(param_path):
            return param_path
            
    # Try to match any params file with the tiff name in it
    param_files = glob.glob(os.path.join(tiff_dir, "params_*.m"))
    for param_file in param_files:
        param_basename = os.path.basename(param_file)
        if tiff_name_no_ext in param_basename:
            return param_file
    
    return None

def random_offset_image_and_update_params(input_path, output_path, canvas_width=155, canvas_height=155):
    """
    Places the input image at a random position within a new black canvas and
    updates the associated parameter file with the new position.
    
    Parameters:
    -----------
    input_path : str
        Path to the input TIFF image
    output_path : str
        Path to save the processed TIFF image
    canvas_width : int
        Width of the new canvas (default: 8000)
    canvas_height : int
        Height of the new canvas (default: 8000)
        
    Returns:
    --------
    tuple or None
        (left_offset_px, top_offset_px) if successful, None otherwise
    """
    try:
        # Open the input image
        img = Image.open(input_path)
        
        # Get original image dimensions
        orig_width, orig_height = img.size
        
        # Calculate maximum possible offsets in pixels
        max_left_offset_px = canvas_width - orig_width
        max_top_offset_px = canvas_height - orig_height
        
        # Check if the image fits within the canvas
        if max_left_offset_px < 0 or max_top_offset_px < 0:
            raise ValueError(f"Original image {input_path} ({orig_width}x{orig_height}) is larger than the canvas ({canvas_width}x{canvas_height})")
        
        # Generate random offsets in pixels
        left_offset_px = random.randint(0, max_left_offset_px)
        top_offset_px = random.randint(0, max_top_offset_px)
        
        # Create a new blank (black) image
        new_img = Image.new(img.mode, (canvas_width, canvas_height), color=0)
        
        # Paste the original image onto the new image at the random position
        new_img.paste(img, (left_offset_px, top_offset_px))
        
        # Save the result
        new_img.save(output_path)
        
        print(f"Processed: {os.path.basename(input_path)} -> {os.path.basename(output_path)}")
        print(f"  Original size: {orig_width}x{orig_height}")
        print(f"  Positioned at ({left_offset_px}, {top_offset_px}) pixels")
        
        return (left_offset_px, top_offset_px)
        
    except Exception as e:
        print(f"Error processing {input_path}: {e}", file=sys.stderr)
        return None

def process_directory(input_dir, output_dir, canvas_width=155, canvas_height=155):
    """
    Processes all TIFF images in the input directory and saves them to the output directory,
    along with updated parameter files.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing TIFF images to process
    output_dir : str
        Directory to save the processed TIFF images and parameter files
    canvas_width : int
        Width of the new canvas (default: 8000)
    canvas_height : int
        Height of the new canvas (default: 8000)
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all TIFF files in the input directory
    tiff_patterns = ['*.tif', '*.tiff', '*.TIF', '*.TIFF']
    tiff_files = []
    for pattern in tiff_patterns:
        tiff_files.extend(glob.glob(os.path.join(input_dir, pattern)))
    
    if not tiff_files:
        print(f"No TIFF files found in {input_dir}")
        return
    
    print(f"Found {len(tiff_files)} TIFF files in {input_dir}")
    
    # Process each TIFF file
    successful_images = 0
    successful_params = 0
    failed_images = 0
    failed_params = 0
    missing_params = 0
    
    for input_path in tiff_files:
        filename = os.path.basename(input_path)
        output_path = os.path.join(output_dir, filename)
        
        # Process the image
        offset_result = random_offset_image_and_update_params(input_path, output_path, canvas_width, canvas_height)
        
        if offset_result:
            successful_images += 1
            left_offset_px, top_offset_px = offset_result
            
            # Find associated parameter file
            param_file_path = find_param_file_for_tiff(input_path)
            
            if param_file_path:
                # Parse the parameter file to get pixel size
                params = parse_parameter_file(param_file_path)
                
                if params and 'pixel_size_nm' in params:
                    # Convert pixel offsets to nanometers
                    pixel_size_nm = params['pixel_size_nm']
                    left_offset_nm = left_offset_px * pixel_size_nm
                    top_offset_nm = top_offset_px * pixel_size_nm
                    
                    # Create output parameter file path
                    param_filename = os.path.basename(param_file_path)
                    output_param_path = os.path.join(output_dir, param_filename)
                    
                    # Update the parameter file
                    if update_parameter_file(param_file_path, output_param_path, left_offset_nm, top_offset_nm):
                        print(f"  Updated parameter file: {param_filename}")
                        print(f"  Position offset in nm: ({left_offset_nm:.2f}, {top_offset_nm:.2f})")
                        successful_params += 1
                    else:
                        print(f"  Failed to update parameter file: {param_filename}")
                        failed_params += 1
                else:
                    # If we can't get pixel size, just copy the parameter file
                    param_filename = os.path.basename(param_file_path)
                    output_param_path = os.path.join(output_dir, param_filename)
                    shutil.copy(param_file_path, output_param_path)
                    print(f"  Warning: Could not extract pixel size from {param_filename}. Parameter file was copied without updates.")
                    failed_params += 1
            else:
                print(f"  No parameter file found for {filename}")
                missing_params += 1
        else:
            failed_images += 1
    
    print(f"\nProcessing complete:")
    print(f"  Images: {successful_images} successful, {failed_images} failed")
    print(f"  Parameter files: {successful_params} updated, {failed_params} failed to update, {missing_params} not found")

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input_dir output_dir", file=sys.stderr)
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Check if input directory exists
    if not os.path.isdir(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    # Process the directory
    process_directory(input_dir, output_dir)

if __name__ == "__main__":
    # Set random seed based on current time for true randomness
    random.seed()
    main()
