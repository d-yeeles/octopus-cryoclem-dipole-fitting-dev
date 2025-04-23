#!/usr/bin/env python3
"""
Script to process all TIFF images in a directory, placing each original image
at a random position within a new 8000×8000 pixel canvas.

Usage:
    python random_offset_tiffs.py input_dir output_dir
    
    - input_dir: directory containing TIFF images to process
    - output_dir: directory to save the processed TIFF images
"""

import os
import sys
import random
from PIL import Image
import glob

def random_offset_image(input_path, output_path, canvas_width=8000, canvas_height=8000):
    """
    Places the input image at a random position within a new black canvas.
    
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
    """
    try:
        # Open the input image
        img = Image.open(input_path)
        
        # Get original image dimensions
        orig_width, orig_height = img.size
        
        # Calculate maximum possible offsets
        max_left_offset = canvas_width - orig_width
        max_top_offset = canvas_height - orig_height
        
        # Check if the image fits within the canvas
        if max_left_offset < 0 or max_top_offset < 0:
            raise ValueError(f"Original image {input_path} ({orig_width}x{orig_height}) is larger than the canvas ({canvas_width}x{canvas_height})")
        
        # Generate random offsets
        left_offset = random.randint(0, max_left_offset)
        top_offset = random.randint(0, max_top_offset)
        
        # Create a new blank (black) image
        new_img = Image.new(img.mode, (canvas_width, canvas_height), color=0)
        
        # Paste the original image onto the new image at the random position
        new_img.paste(img, (left_offset, top_offset))
        
        # Save the result
        new_img.save(output_path)
        
        print(f"Processed: {os.path.basename(input_path)} -> {os.path.basename(output_path)}")
        print(f"  Original size: {orig_width}x{orig_height}")
        print(f"  Positioned at ({left_offset}, {top_offset})")
        
        return True
        
    except Exception as e:
        print(f"Error processing {input_path}: {e}", file=sys.stderr)
        return False

def process_directory(input_dir, output_dir, canvas_width=8000, canvas_height=8000):
    """
    Processes all TIFF images in the input directory and saves them to the output directory.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing TIFF images to process
    output_dir : str
        Directory to save the processed TIFF images
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
    successful = 0
    failed = 0
    for input_path in tiff_files:
        filename = os.path.basename(input_path)
        output_path = os.path.join(output_dir, filename)
        
        if random_offset_image(input_path, output_path, canvas_width, canvas_height):
            successful += 1
        else:
            failed += 1
    
    print(f"\nProcessing complete: {successful} successful, {failed} failed")

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
