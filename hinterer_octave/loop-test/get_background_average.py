import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Load the image as a NumPy array
# Replace 'your_image.jpg' with your image file path
image_path = 'ASIL240923C05_50.ome.tif-frames/middle_frame_01.tif'
image = np.array(Image.open(image_path).convert('L'))  # Convert to grayscale

# Define the patch coordinates (row_start:row_end, col_start:col_end)
row_start, row_end = 75, 125  # Rows 100 to 150
col_start, col_end = 75, 125  # Columns 200 to 250

# Calculate the mean value of the patch
patch = image[row_start:row_end, col_start:col_end]
patch_mean = np.mean(patch)
print(f"Mean value of the patch: {patch_mean}")

# Plot the image with the patch outlined
fig, ax = plt.subplots()
ax.imshow(image, cmap='gray')
ax.set_title('Image with Patch Highlighted')

# Add a red rectangle to indicate the patch
rect = patches.Rectangle(
    (col_start, row_start),  # Bottom-left corner
    col_end - col_start,     # Width
    row_end - row_start,     # Height
    linewidth=2, edgecolor='red', facecolor='none'
)
ax.add_patch(rect)

# Show the plot
plt.show()

