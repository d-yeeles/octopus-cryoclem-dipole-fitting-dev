from PIL import Image
import numpy as np

# Load the TIF stack
tif_stack = Image.open('stack_all_simulations.tif')

# Create a list to store rotated frames
rotated_frames = []

# Process each frame in the stack
try:
    frame_count = 0
    while True:
        # Select the frame
        tif_stack.seek(frame_count)
        
        # Convert to numpy array, rotate, and convert back to PIL Image
        frame = np.array(tif_stack)
        rotated_frame = np.rot90(frame, 2)  # Rotate 180 degrees
        rotated_frames.append(Image.fromarray(rotated_frame))
        
        frame_count += 1
except EOFError:
    # End of frames
    pass

# Save the first frame
rotated_frames[0].save(
    'stack_all_simulations_rotated.tif',
    save_all=True,
    append_images=rotated_frames[1:],  # Append the rest of the frames
    compression='tiff_deflate',
    compression_level=9
)

print(f"Rotated {frame_count} frames and saved to stack_all_simulations_rotated.tif")
