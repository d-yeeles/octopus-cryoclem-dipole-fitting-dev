#!/usr/bin/env python3
"""
Script to add empty (black) pixels around a TIFF image, 
positioning the original image off-center.

Usage:
    python pad_tiff.py input.tif output.tif left_pad top_pad new_width new_height
    
    - input.tif: path to the input TIFF image
    - output.tif: path to save the padded TIFF image
    - left_pad: number of pixels to pad on the left side
    - top_pad: number of pixels to pad on the top
    - new_width: total width of the new image (must be >= original width + left_pad)
    - new_height: total height of the new image (must be >= original height + top_pad)
"""

import sys
from PIL import Image
import os

def pad_image_off_center(input_path, output_path, left_pad, top_pad, new_width, new_height):
    """
    Adds empty (black) pixels around a TIFF image, positioning the original image off-center.
    
    Parameters:
    -----------
    input_path : str
        Path to the input TIFF image
    output_path : str
        Path to save the padded TIFF image
    left_pad : int
        Number of pixels to pad on the left side
    top_pad : int
        Number of pixels to pad on the top
    new_width : int
        Total width of the new image
    new_height : int
        Total height of the new image
    """
    try:
        # Open the input image
        img = Image.open(input_path)
        
        # Get original image dimensions
        orig_width, orig_height = img.size
        
        # Check if the new dimensions are sufficient
        if new_width < orig_width + left_pad:
            raise ValueError(f"New width ({new_width}) must be at least original width + left padding ({orig_width + left_pad})")
        if new_height < orig_height + top_pad:
            raise ValueError(f"New height ({new_height}) must be at least original height + top padding ({orig_height + top_pad})")
        
        # Create a new blank (black) image
        new_img = Image.new(img.mode, (new_width, new_height), color=0)
        
        # Paste the original image onto the new image at the specified position
        new_img.paste(img, (left_pad, top_pad))
        
        # Save the result
        new_img.save(output_path)
        
        print(f"Successfully created padded image: {output_path}")
        print(f"Original size: {orig_width}x{orig_height}")
        print(f"New size: {new_width}x{new_height}")
        print(f"Image positioned at ({left_pad}, {top_pad})")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    if len(sys.argv) != 7:
        print(f"Usage: {sys.argv[0]} input.tif output.tif left_pad top_pad new_width new_height", file=sys.stderr)
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    left_pad = int(sys.argv[3])
    top_pad = int(sys.argv[4])
    new_width = int(sys.argv[5])
    new_height = int(sys.argv[6])
    
    pad_image_off_center(input_path, output_path, left_pad, top_pad, new_width, new_height)

if __name__ == "__main__":
    main()
