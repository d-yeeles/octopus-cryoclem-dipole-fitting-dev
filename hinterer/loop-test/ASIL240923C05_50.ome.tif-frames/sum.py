from PIL import Image
import numpy as np

image = Image.open("middle_frame_02.tif")
data = np.array(image)
photon_count = np.sum(data)

print("Total Photon Count:", photon_count)

