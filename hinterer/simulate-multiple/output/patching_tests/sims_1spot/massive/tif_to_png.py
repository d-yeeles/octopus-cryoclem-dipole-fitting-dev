#!/usr/bin/env python3

import os
import sys
import numpy as np
from PIL import Image
import tifffile
import glob
from tqdm import tqdm

def normalize_image(image):
    """
    Normalize image exactly like the MATLAB code:
    (image - min) / (max - min)
    """
    # Convert to float64 for calculations
    img_float = image.astype(np.float64)
    
    # Find min and max values
    min_val = np.min(img_float)
    max_val = np.max(img_float)
    
    # Avoid division by zero
    if max_val == min_val:
        return np.zeros_like(img_float)
    
    # Normalize to 0-1 range exactly like the MATLAB code
    normalized = (img_float - min_val) / (max_val - min_val)
    
    # Convert to 8-bit (0-255)
    img_8bit = (normalized * 255).astype(np.uint8)
    
    return img_8bit

def convert_tif_to_png(input_dir, output_dir):
    """Convert all TIFF files in input_dir to PNG in output_dir."""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get all TIFF files
    tif_files = glob.glob(os.path.join(input_dir, "*.tif")) + glob.glob(os.path.join(input_dir, "*.tiff"))
    
    print(f"Found {len(tif_files)} TIFF files to convert")
    
    # Process each file
    for tif_file in tqdm(tif_files, desc="Converting"):
        try:
            # Get base filename without extension
            base_name = os.path.splitext(os.path.basename(tif_file))[0]
            output_file = os.path.join(output_dir, f"{base_name}.png")
            
            # Read the TIFF file with tifffile which handles 32-bit TIFFs well
            image = tifffile.imread(tif_file)
            
            # Print image info for debugging
            print(f"\nFile: {tif_file}")
            print(f"Shape: {image.shape}, Dtype: {image.dtype}")
            print(f"Min: {np.min(image)}, Max: {np.max(image)}")
            
            # Normalize the image
            normalized = normalize_image(image)
            
            # Check if normalized image is all black (or nearly all black)
            black_percentage = np.sum(normalized < 5) / normalized.size * 100
            if black_percentage > 95:
                print(f"  Warning: Image is {black_percentage:.1f}% black. Trying alternative normalization...")
                
                # Try alternative normalization method with log scaling
                # This can help if the data has a few extreme values
                eps = 1e-10  # Small value to avoid log(0)
                log_image = np.log(image.astype(np.float64) + eps)
                normalized = normalize_image(log_image)
                
                # Check again
                black_percentage = np.sum(normalized < 5) / normalized.size * 100
                if black_percentage > 99.999:
                    print(f"  Warning: Image is still {black_percentage:.1f}% black after log normalization.")
                    
                    # Try with histogram equalization
                    from skimage import exposure
                    p2, p98 = np.percentile(image, (2, 98))
                    img_rescale = exposure.rescale_intensity(image, in_range=(p2, p98))
                    normalized = normalize_image(img_rescale)
                    
                    black_percentage = np.sum(normalized < 5) / normalized.size * 100
                    print(f"  After histogram equalization: {black_percentage:.1f}% black")
            
            # Save as PNG
            Image.fromarray(normalized).save(output_file)
            print(f"  Saved to {output_file}")
            
        except Exception as e:
            print(f"Error processing {tif_file}: {e}")

if __name__ == "__main__":
    # Get input and output directories from arguments or use defaults
    input_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    output_dir = os.path.join(input_dir, "png_converted")
    
    convert_tif_to_png(input_dir, output_dir)
    print(f"\nConversion complete. PNG files are in {output_dir}")
