#!/usr/bin/env python3
import numpy as np
from PIL import Image
import os

def pad_tif_with_zeros(input_path, output_path=None, padding=100):
    """
    Pad a TIF image with zeros (black) on all sides.
    
    Parameters:
    -----------
    input_path : str
        Path to the input TIF file
    output_path : str, optional
        Path to save the padded TIF file. If None, will use '[original_name]_padded.tif'
    padding : int, optional
        Number of pixels to pad on each side (default: 100)
    
    Returns:
    --------
    str
        Path to the saved padded image
    """
    # Open the TIF file
    img = Image.open(input_path)
    img_array = np.array(img)
    
    # Get original dimensions
    if len(img_array.shape) == 2:
        # Grayscale image
        height, width = img_array.shape
        channels = 1
    else:
        # Color image
        height, width, channels = img_array.shape
    
    # Create new array with padding
    if channels == 1:
        padded_array = np.zeros((height + 2 * padding, width + 2 * padding), dtype=img_array.dtype)
        # Place the original image in the center
        padded_array[padding:padding + height, padding:padding + width] = img_array
    else:
        padded_array = np.zeros((height + 2 * padding, width + 2 * padding, channels), dtype=img_array.dtype)
        # Place the original image in the center
        padded_array[padding:padding + height, padding:padding + width, :] = img_array
    
    # Create a new image from the padded array
    padded_img = Image.fromarray(padded_array)
    
    # If output path is not provided, create one
    if output_path is None:
        base, ext = os.path.splitext(input_path)
        output_path = f"{base}_padded{ext}"
    
    # Save the padded image
    padded_img.save(output_path)
    
    print(f"Image padded and saved to {output_path}")
    return output_path

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Pad a TIF image with zeros (black) on all sides")
    parser.add_argument("input_file", help="Path to the input TIF file")
    parser.add_argument("-o", "--output", help="Path to save the padded TIF file (optional)")
    parser.add_argument("-p", "--padding", type=int, default=100, help="Number of pixels to pad on each side (default: 100)")
    
    args = parser.parse_args()
    
    pad_tif_with_zeros(args.input_file, args.output, args.padding)
