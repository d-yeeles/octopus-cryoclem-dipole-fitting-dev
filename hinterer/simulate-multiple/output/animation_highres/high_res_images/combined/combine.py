#!/usr/bin/env python3
import os
from PIL import Image, ImageDraw, ImageFont

def combine_images(left_img_path, right_img_path, theta_number):
    # Open the images
    left_img = Image.open(left_img_path)
    right_img = Image.open(right_img_path)
    
    # Ensure images have the same height
    max_height = max(left_img.height, right_img.height)
    
    # Resize images to have the same height while maintaining aspect ratio
    left_img = resize_image(left_img, max_height)
    right_img = resize_image(right_img, max_height)
    
    # Calculate total width
    total_width = left_img.width + right_img.width
    
    # Create a new blank image
    combined_img = Image.new('RGB', (total_width, max_height), color='white')
    
    # Paste the left and right images
    combined_img.paste(left_img, (0, 0))
    combined_img.paste(right_img, (left_img.width, 0))
    
    # Add theta number text
    draw = ImageDraw.Draw(combined_img)
    
    # Try to use a default font, with fallback
    try:
        # Attempt to use a standard font (now with 1/4 the original size)
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 48)
    except IOError:
        # Fallback to default font
        font = ImageFont.load_default()
    
    # Remove leading zero and format text
    theta_display = str(int(theta_number))
    text = f"θ = {theta_display}°"
    
    # Use getbbox() for text size in newer Pillow versions
    left, top, right, bottom = draw.textbbox((0, 0), text, font=font)
    text_width = right - left
    text_height = bottom - top
    
    # Position for text (bottom center)
    text_position = ((total_width - text_width) // 2, max_height - text_height - 15)
    
    # Draw text with white color and black outline for visibility
    # Draw outline
    outline_color = (0, 0, 0)  # Black outline
    for offset in [(1,1), (-1,1), (1,-1), (-1,-1)]:
        draw.text((text_position[0]+offset[0], text_position[1]+offset[1]), 
                  text, font=font, fill=outline_color)
    
    # Draw white text on top
    draw.text(text_position, text, fill=(255,255,255), font=font)
    
    # Save the combined image
    output_filename = f"theta{theta_number}.png"
    combined_img.save(output_filename)
    print(f"Created {output_filename}")
    
    # Close the images
    left_img.close()
    right_img.close()
    combined_img.close()

def resize_image(img, target_height):
    # Calculate the new width while maintaining aspect ratio
    aspect_ratio = img.width / img.height
    new_width = int(target_height * aspect_ratio)
    return img.resize((new_width, target_height), Image.LANCZOS)

def main():
    # Create dictionaries to store image paths
    hinterer_images = {}
    mortensen_images = {}
    
    # Collect image paths
    for filename in os.listdir('.'):
        if filename.startswith('hinterer_theta') and filename.endswith('.png'):
            theta = filename.split('theta')[1].split('.png')[0]
            hinterer_images[theta] = filename
        elif filename.startswith('mortensen_theta') and filename.endswith('.png'):
            theta = filename.split('theta')[1].split('.png')[0]
            mortensen_images[theta] = filename
    
    # Combine matching images
    for theta in sorted(hinterer_images.keys()):
        if theta in mortensen_images:
            left_path = hinterer_images[theta]
            right_path = mortensen_images[theta]
            combine_images(left_path, right_path, theta)

if __name__ == '__main__':
    main()
