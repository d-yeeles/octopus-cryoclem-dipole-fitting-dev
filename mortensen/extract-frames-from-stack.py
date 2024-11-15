import tifffile
import numpy as np
import os

# Ask the user for input
input_file = input("Enter the name of the TIFF file (e.g., ome.tif): ")
num_frames = int(input("Enter the number of frames to extract (e.g., 10): "))

num_frames = num_frames*2

# Read the TIFF stack using tifffile
with tifffile.TiffFile(input_file) as tif:
    total_frames = len(tif.pages)
    print(f"Total frames in the TIFF file: {total_frames}")

    # Check if the number of frames requested is less than or equal to the total frames
    if total_frames < num_frames:
        print(f"The TIFF file contains only {total_frames} frames, which is less than the number of frames you want to extract.")
        exit(1)

    # Calculate the middle frames
    middle_start = total_frames // 2 - num_frames // 2
    middle_end = total_frames // 2 + num_frames // 2

    # Ensure the end index doesn't exceed the total number of frames
    if middle_end > total_frames:
        middle_end = total_frames

    # Display the frame range
    print(f"Extracting frames from {middle_start} to {middle_end - 1}.")

    # Create an output folder to save the frames
    output_dir = input_file + "-frames"
    os.makedirs(output_dir, exist_ok=True)

    # Extract every other frame from the middle frames
    extracted_frame_count = 0
    for i in range(middle_start, middle_end, 2):  # Skip every other frame
        # Read the frame
        frame = tif.pages[i].asarray()

        # Save the frame to a TIFF file
        output_filename = os.path.join(output_dir, f"middle_frame_{extracted_frame_count + 1:02d}.tif")
        tifffile.imwrite(output_filename, frame)

        # Increment frame counter
        extracted_frame_count += 1

    print(f"Extracted {extracted_frame_count} frames to {output_dir}.")
