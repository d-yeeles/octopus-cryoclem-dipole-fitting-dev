from PIL import Image, ImageDraw, ImageFont
import os

# Directory containing your images
img_dir = "./mortensen"  # Change if the images are in a different folder

# Font settings
font_size = 50
font_color = "white"
font_path = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"  # A common font on Ubuntu

# Get sorted file list
files = sorted([f for f in os.listdir(img_dir) if f.startswith("sim_frame") and f.endswith(".png")])

for i, file_name in enumerate(files):
    number = i * 10  # 0, 10, 20, ...
    img_path = os.path.join(img_dir, file_name)
    
    with Image.open(img_path) as img:
        draw = ImageDraw.Draw(img)
        try:
            font = ImageFont.truetype(font_path, font_size)
        except IOError:
            font = ImageFont.load_default()
        
        text = str(number)+'°'
        bbox = draw.textbbox((0, 0), text, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]
        
        # Padding from the edge
        padding = 10
        x = img.width - text_width - padding
        y = img.height - text_height - padding
        
        draw.text((x, y), text, font=font, fill=font_color)
        img.save(os.path.join(img_dir, f"numbered_{file_name}"))

print("Done! Numbered images saved with 'numbered_' prefix.")

