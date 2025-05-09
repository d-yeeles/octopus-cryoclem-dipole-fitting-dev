import os
import numpy as np
from tifffile import imread, imwrite

def stack_tif_stacks(stack_filenames, output_filename="combined_stack.tif"):
    """
    Combine multiple TIF stacks into a single stack.
    
    Parameters:
    -----------
    stack_filenames : list
        List of filenames of TIF stacks to combine
    output_filename : str, optional
        Name of the output combined stack file
    """
    print(f"Processing {len(stack_filenames)} TIF stacks...")
    
    # List to hold all images from all stacks
    all_images = []
    failed_stacks = []
    
    # Load each stack and append its images to all_images
    for i, stack_file in enumerate(stack_filenames):
        try:
            # Load the stack
            print(f"Loading stack {i+1}/{len(stack_filenames)}: {stack_file}")
            stack = imread(stack_file)
            
            # If the stack is a single image, reshape it to have a dimension for frames
            if len(stack.shape) == 2:
                stack = np.expand_dims(stack, axis=0)
            
            # Get number of frames in this stack
            num_frames = stack.shape[0]
            print(f"  - Found {num_frames} frames in {stack_file}")
            
            # Append each frame to all_images
            for j in range(num_frames):
                all_images.append(stack[j])
                
            print(f"  - Successfully added all frames from {stack_file}")
            
        except Exception as e:
            print(f"Error loading {stack_file}: {e}")
            failed_stacks.append(stack_file)
    
    # Report summary
    print(f"Successfully loaded {len(all_images)} frames from {len(stack_filenames) - len(failed_stacks)} stacks.")
    if failed_stacks:
        print(f"Failed to load {len(failed_stacks)} stacks: {', '.join(failed_stacks)}")
    
    # If we have images, save as a multi-page TIFF
    if all_images:
        # Convert list of images to a 3D numpy array
        combined_stack = np.array(all_images)
        
        # Save as multi-page TIFF
        try:
            imwrite(output_filename, combined_stack, compression='zlib', photometric='minisblack')
            print(f"Successfully created combined TIFF stack with {len(all_images)} frames at {output_filename}")
        except Exception as e:
            print(f"Error saving combined stack: {e}")
    else:
        print("No images were successfully loaded. Check file paths and image format.")

if __name__ == "__main__":

    stack_files = [f for f in os.listdir('.') if f.endswith('.tif') and os.path.isfile(f)]
    print(f"Found {len(stack_files)} .tif files in current directory")
    
    # Set the output filename
    output_file = "combined_sim_stack.tif"
    
    # Run the stacking function
    stack_tif_stacks(stack_files, output_file)
