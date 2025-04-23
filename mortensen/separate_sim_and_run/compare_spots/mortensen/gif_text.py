import os
from PIL import Image, ImageDraw, ImageFont
import glob

def add_text_to_images(input_dir, output_dir, font_size=36):
    """Add sequential numbering text to PNG images and save to output directory"""
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Get all PNG files in the directory
    png_files = sorted(glob.glob(os.path.join(input_dir, "*.png")))
    
    # Check if we have enough files
    if len(png_files) < 91:
        print(f"Warning: Only found {len(png_files)} PNG files. Will number them sequentially.")
    
    # Ubuntu 22.04 common font paths
    font_paths = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/ubuntu/Ubuntu-R.ttf",
        "/usr/share/fonts/truetype/freefont/FreeSans.ttf",
        "/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf",
        # Fallbacks
        "DejaVuSans.ttf",
        "arial.ttf"
    ]
    
    font = None
    for font_path in font_paths:
        try:
            font = ImageFont.truetype(font_path, font_size)
            print(f"Using font: {font_path} with size {font_size}")
            break
        except IOError:
            continue
    
    # If no font was loaded, use the default font
    if font is None:
        print("Could not load any TrueType fonts. Using default font.")
        font = ImageFont.load_default()
    
    # Process each image
    for i, file_path in enumerate(png_files[:91]):  # Limit to 91 files (0-90)
        img = Image.open(file_path)
        
        # Convert the image to RGB mode if it's not already
        if img.mode != 'RGB':
            img = img.convert('RGB')
        
        # Create a drawing object
        draw = ImageDraw.Draw(img)
        
        # Add text - white with slight shadow for better visibility
        text = str(i)
        
        # Get text dimensions - use a simple approximation if methods fail
        try:
            left, top, right, bottom = draw.textbbox((0, 0), text, font=font)
            text_width = right - left
            text_height = bottom - top
        except (AttributeError, TypeError):
            # Simple approximation based on font size
            text_width = len(text) * font_size * 0.6
            text_height = font_size * 1.2
        
        # Position text in bottom right with some padding
        position = (img.width - int(text_width) - 20, img.height - int(text_height) - 20)
        
        # Colors for text
        black = (0, 0, 0)
        white = (255, 255, 255)
        
        # Add shadow (black outline) for better visibility
        shadow_offset = max(1, font_size // 20)  # Scale shadow with font size
        shadow_positions = [
            (position[0]-shadow_offset, position[1]-shadow_offset), 
            (position[0]-shadow_offset, position[1]+shadow_offset),
            (position[0]+shadow_offset, position[1]-shadow_offset),
            (position[0]+shadow_offset, position[1]+shadow_offset)
        ]
        
        for shadow_pos in shadow_positions:
            draw.text(shadow_pos, text, fill=black, font=font)
        
        # Draw the white text
        draw.text(position, text, fill=white, font=font)
        
        # Save the modified image
        output_path = os.path.join(output_dir, f"numbered_{i:03d}.png")
        img.save(output_path)
        print(f"Processed image {i}: {output_path}")
    
    return len(png_files[:91])

def create_animated_gif(input_dir, output_file, duration=100):
    """Create an animated GIF from numbered PNG files with bounce effect"""
    # Get all numbered PNG files
    numbered_files = sorted(glob.glob(os.path.join(input_dir, "numbered_*.png")))
    
    if not numbered_files:
        print("No numbered PNG files found!")
        return False
    
    # Open all images
    images = [Image.open(file) for file in numbered_files]
    
    # Create the forward sequence
    forward_frames = images.copy()
    
    # Create the reverse sequence (excluding first and last frame to avoid duplicates)
    reverse_frames = images[-2:0:-1]  # All frames except first and last, in reverse order
    
    # Combine forward and reverse sequences for bounce effect
    all_frames = forward_frames + reverse_frames
    
    # Save as animated GIF
    all_frames[0].save(
        output_file,
        format="GIF",
        append_images=all_frames[1:],
        save_all=True,
        duration=duration,  # Duration per frame in milliseconds
        loop=0  # 0 means loop indefinitely
    )
    
    print(f"Created animated GIF: {output_file}")
    return True

def main():
    # Directories and output file
    input_dir = "input_pngs"  # Directory containing the original PNG files
    output_dir = "numbered_pngs"  # Directory to save numbered PNG files
    output_gif = "animated_sequence.gif"  # Output GIF file
    
    # Get user input for directories
    user_input_dir = input(f"Enter input PNG directory (default: {input_dir}): ")
    if user_input_dir:
        input_dir = user_input_dir
    
    user_output_dir = input(f"Enter output directory for numbered PNGs (default: {output_dir}): ")
    if user_output_dir:
        output_dir = user_output_dir
    
    user_output_gif = input(f"Enter output GIF filename (default: {output_gif}): ")
    if user_output_gif:
        output_gif = user_output_gif
    
    # Get font size from user
    font_size = input("Enter font size (default: 36): ")
    try:
        font_size = int(font_size) if font_size else 36
    except ValueError:
        font_size = 36
        print("Invalid font size. Using default: 36")
    
    # Add numbering to images
    num_processed = add_text_to_images(input_dir, output_dir, font_size)
    print(f"Processed {num_processed} images with numbering")
    
    # Create animated GIF
    gif_duration = int(input("Enter frame duration in milliseconds (default: 100): ") or "100")
    success = create_animated_gif(output_dir, output_gif, gif_duration)
    
    if success:
        print(f"Process complete! Animated GIF saved as {output_gif}")
    else:
        print("Failed to create animated GIF")

if __name__ == "__main__":
    main()
