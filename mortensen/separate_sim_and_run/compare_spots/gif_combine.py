from PIL import Image
import os
import sys

def combine_gifs(gif1_path, gif2_path, output_path):
    """
    Combine two animated GIFs side by side into a new animated GIF.
    
    Parameters:
    gif1_path (str): Path to the first GIF file
    gif2_path (str): Path to the second GIF file
    output_path (str): Path where the combined GIF will be saved
    """
    # Open the GIFs
    gif1 = Image.open(gif1_path)
    gif2 = Image.open(gif2_path)
    
    # Get the dimensions
    width1, height1 = gif1.size
    width2, height2 = gif2.size
    
    # Calculate the dimensions of the combined GIF
    new_width = width1 + width2
    new_height = max(height1, height2)
    
    # Extract all frames from both GIFs
    frames1 = []
    frames2 = []
    
    try:
        # Extract frames from first GIF
        while True:
            # Get the current frame and its duration
            current = gif1.copy()
            frames1.append({
                'image': current,
                'duration': gif1.info.get('duration', 100)
            })
            gif1.seek(gif1.tell() + 1)
    except EOFError:
        pass  # End of frames for the first GIF
    
    try:
        # Extract frames from second GIF
        while True:
            # Get the current frame and its duration
            current = gif2.copy()
            frames2.append({
                'image': current,
                'duration': gif2.info.get('duration', 100)
            })
            gif2.seek(gif2.tell() + 1)
    except EOFError:
        pass  # End of frames for the second GIF
    
    # Determine the number of frames in the combined GIF
    # (we'll loop the shorter GIF if necessary)
    num_frames = max(len(frames1), len(frames2))
    
    # Create the combined frames
    combined_frames = []
    durations = []
    
    for i in range(num_frames):
        # Get frames from both GIFs (loop if necessary)
        frame1_idx = i % len(frames1)
        frame2_idx = i % len(frames2)
        
        frame1 = frames1[frame1_idx]['image']
        frame2 = frames2[frame2_idx]['image']
        
        # Choose the shorter duration to keep animations in sync
        duration = min(frames1[frame1_idx]['duration'], frames2[frame2_idx]['duration'])
        durations.append(duration)
        
        # Create a new blank frame
        new_frame = Image.new('RGBA', (new_width, new_height))
        
        # Paste both frames side by side
        new_frame.paste(frame1, (0, 0))
        new_frame.paste(frame2, (width1, 0))
        
        # Convert to 'P' mode (palette) for GIF compatibility
        if new_frame.mode != 'P':
            new_frame = new_frame.convert('P')
        
        combined_frames.append(new_frame)
    
    # Save the combined GIF
    if combined_frames:
        combined_frames[0].save(
            output_path,
            save_all=True,
            append_images=combined_frames[1:],
            optimize=False,
            duration=durations,
            loop=0  # 0 means loop indefinitely
        )
        print(f"Combined GIF saved to {output_path}")
    else:
        print("No frames were extracted from the input GIFs.")

if __name__ == "__main__":
    # Check if command line arguments are provided
    if len(sys.argv) == 4:
        gif1_path = sys.argv[1]
        gif2_path = sys.argv[2]
        output_path = sys.argv[3]
    else:
        # Otherwise, prompt for input
        gif1_path = input("Enter the path to the first GIF: ")
        gif2_path = input("Enter the path to the second GIF: ")
        output_path = input("Enter the output path for the combined GIF: ")
    
    # Verify files exist
    if not os.path.exists(gif1_path):
        print(f"Error: File {gif1_path} does not exist.")
        sys.exit(1)
    if not os.path.exists(gif2_path):
        print(f"Error: File {gif2_path} does not exist.")
        sys.exit(1)
    
    # Combine the GIFs
    combine_gifs(gif1_path, gif2_path, output_path)
