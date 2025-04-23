import os
import argparse
from tifffile import imread, imwrite
import numpy as np

def flip_tif_stack(input_path, output_path=None):
    """
    Flip a TIF stack horizontally.
    
    Parameters:
    -----------
    input_path : str
        Path to the input TIF stack file
    output_path : str, optional
        Path to save the flipped TIF stack file. If None, a new filename will be generated.
    """
    # Generate output path if not provided
    if output_path is None:
        base, ext = os.path.splitext(input_path)
        output_path = f"{base}_flipped{ext}"
    
    # Read the TIF stack
    print(f"Reading TIF stack from {input_path}...")
    stack = imread(input_path)
    
    # Check dimensions and print info
    print(f"Stack shape: {stack.shape}")
    
    # Flip the stack horizontally
    # For a 3D stack, axis=2 flips each frame horizontally
    # For a 2D image, axis=1 flips horizontally
    if len(stack.shape) == 3:
        flipped_stack = np.flip(stack, axis=2)
        print(f"Flipped 3D stack along horizontal axis")
    else:
        flipped_stack = np.flip(stack, axis=1)
        print(f"Flipped 2D image along horizontal axis")
        
    # Save the flipped stack
    print(f"Saving flipped stack to {output_path}...")
    imwrite(output_path, flipped_stack)
    print("Done!")
    
    return output_path

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Flip a TIF stack horizontally")
    parser.add_argument("input", help="Path to input TIF stack file")
    parser.add_argument("-o", "--output", help="Path to output flipped TIF stack file (optional)")
    args = parser.parse_args()
    
    # Flip the stack
    output_file = flip_tif_stack(args.input, args.output)
    print(f"Successfully flipped TIF stack and saved to {output_file}")
