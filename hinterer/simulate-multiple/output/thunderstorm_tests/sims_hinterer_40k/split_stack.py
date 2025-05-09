import os
import numpy as np
from tifffile import imread, imwrite

def split_tif_stack(input_filename, num_stacks=10, output_prefix="split_stack_"):
    """
    Split a single TIF stack into multiple smaller stacks.
    
    Parameters:
    -----------
    input_filename : str
        Filename of the TIF stack to split
    num_stacks : int, optional
        Number of stacks to split into
    output_prefix : str, optional
        Prefix for output filenames
    """
    print(f"Processing {input_filename} to create {num_stacks} TIF stacks...")
    
    try:
        # Load the input stack
        print(f"Loading stack: {input_filename}")
        stack = imread(input_filename)
        
        # If the stack is a single image, reshape it to have a dimension for frames
        if len(stack.shape) == 2:
            stack = np.expand_dims(stack, axis=0)
        
        # Get number of frames in this stack
        total_frames = stack.shape[0]
        print(f"  - Found {total_frames} frames in {input_filename}")
        
        # Calculate frames per output stack
        frames_per_stack = total_frames // num_stacks
        remainder = total_frames % num_stacks
        
        print(f"  - Will create {num_stacks} stacks with approximately {frames_per_stack} frames each")
        
        # Create output stacks
        start_idx = 0
        for i in range(num_stacks):
            # Calculate how many frames to put in this stack
            # Distribute remainder frames across first 'remainder' stacks
            current_frames = frames_per_stack + (1 if i < remainder else 0)
            end_idx = start_idx + current_frames
            
            # Check if we're at the end of the input stack
            if start_idx >= total_frames:
                break
                
            # Create output filename
            output_filename = f"{output_prefix}{i+1:02d}.tif"
            
            # Extract frames for this stack
            output_stack = stack[start_idx:end_idx]
            
            # Save as multi-page TIFF
            try:
                imwrite(output_filename, output_stack, compression='zlib', photometric='minisblack')
                print(f"  - Created stack {i+1}/{num_stacks}: {output_filename} with {current_frames} frames")
            except Exception as e:
                print(f"  - Error saving {output_filename}: {e}")
            
            # Update start index for next stack
            start_idx = end_idx
        
        print(f"Successfully split {input_filename} into {num_stacks} stacks")
        
    except Exception as e:
        print(f"Error processing {input_filename}: {e}")

if __name__ == "__main__":
    # Specify your input stack file
    input_file = "sim_stack.tif"
    
    # Set the number of output stacks
    num_output_stacks = 10
    
    # Set the output filename prefix
    output_prefix = "split_stack_"
    
    # Run the splitting function
    split_tif_stack(input_file, num_output_stacks, output_prefix)
