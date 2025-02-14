import argparse
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

def view_tiff_image(file_path):
    # Open the TIFF image file
    try:
        tiff_image = Image.open(file_path)
        # Convert image to numpy array
        image_data = np.array(tiff_image, dtype=np.float64)  # Convert to double precision
        
        # Normalize the image (equivalent to MATLAB's normalization)
        image_data = (image_data - np.min(image_data)) / (np.max(image_data) - np.min(image_data))
        
        # Display the image
        plt.imshow(image_data, cmap='gray')  # Display as grayscale
        plt.axis('off')  # Hide the axis
        plt.show()
    except Exception as e:
        print(f"Error opening image: {e}")

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="View a 32-bit TIFF image from a file path.")
    parser.add_argument("file_path", type=str, help="Path to the TIFF image file")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # View the image using the provided file path
    view_tiff_image(args.file_path)

