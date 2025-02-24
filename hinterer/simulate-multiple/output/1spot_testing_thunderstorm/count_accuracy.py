import pandas as pd
from PIL import Image

# Load the CSV file
csv_file = "thunderstorm_results.csv"  # Change this to your actual file path
df = pd.read_csv(csv_file)

# Determine the number of frames in the TIFF stack
tiff_file = "all_stack.tif"  # Change this to your actual TIFF file
with Image.open(tiff_file) as img:
    num_frames = img.n_frames

# Define polar and azimuth angle steps
polar_angles = list(range(0, 91, 1))
azimuth_angles = list(range(0, 361, 1))

# Generate frame-to-angle mapping
frame_to_angle = {}
frame_index = 1
for polar in polar_angles:
    for azimuth in azimuth_angles:
        frame_to_angle[frame_index] = (polar, azimuth)
        frame_index += 1

# Expected range of frames
expected_frames = set(frame_to_angle.keys())

# Get actual frames from the CSV
actual_frames = df["frame"].tolist()

# Identify missing and duplicate frames
missing_frames = sorted(expected_frames - set(actual_frames))
duplicate_frames = sorted(set(f for f in actual_frames if actual_frames.count(f) > 1))
correct_frames = sorted(expected_frames - set(missing_frames) - set(duplicate_frames))

# Extract angles for each category
missing_inc = [frame_to_angle[f][0] for f in missing_frames]
missing_az = [frame_to_angle[f][1] for f in missing_frames]
duplicate_inc = [frame_to_angle[f][0] for f in duplicate_frames]
duplicate_az = [frame_to_angle[f][1] for f in duplicate_frames]
correct_inc = [frame_to_angle[f][0] for f in correct_frames]
correct_az = [frame_to_angle[f][1] for f in correct_frames]

# Write results to a file
output_file = "accuracy_results.py"
with open(output_file, "w") as f:
    f.write(f"false_negatives = {missing_frames}\n")
    f.write(f"false_negatives_inc = {missing_inc}\n")
    f.write(f"false_negatives_az = {missing_az}\n")
    f.write(f"false_positives = {duplicate_frames}\n")
    f.write(f"false_positives_inc = {duplicate_inc}\n")
    f.write(f"false_positives_az = {duplicate_az}\n")
    f.write(f"correct = {correct_frames}\n")
    f.write(f"correct_inc = {correct_inc}\n")
    f.write(f"correct_az = {correct_az}\n")

print(f"Results saved to {output_file}")


# -----------------------------------------


import matplotlib.pyplot as plt
import numpy as np
from accuracy_results import (
    false_negatives_az, false_negatives_inc, 
    false_positives_az, false_positives_inc, 
    correct_az, correct_inc
)

# Define bin edges based on azimuth and inclination steps
azimuth_bins = np.arange(0, 361, 4)
inclination_bins = np.arange(0, 91, 1)

# Create subplots
fig, axes = plt.subplots(3, 2, figsize=(12, 12))

# Plot false positives
axes[0, 0].hist(false_positives_az, bins=azimuth_bins, edgecolor='black', alpha=0.7)
axes[0, 0].set_title("False Positives")
axes[0, 0].set_xlabel("φ, °")
axes[0, 0].set_ylabel("Count")

axes[0, 1].hist(false_positives_inc, bins=inclination_bins, edgecolor='black', alpha=0.7)
axes[0, 1].set_title("False Positives")
axes[0, 1].set_xlabel("θ, °")
axes[0, 1].set_ylabel("Count")

# Plot false negatives
axes[1, 0].hist(false_negatives_az, bins=azimuth_bins, edgecolor='black', alpha=0.7)
axes[1, 0].set_title("False Negatives")
axes[1, 0].set_xlabel("φ, °")
axes[1, 0].set_ylabel("Count")

axes[1, 1].hist(false_negatives_inc, bins=inclination_bins, edgecolor='black', alpha=0.7)
axes[1, 1].set_title("False Negatives")
axes[1, 1].set_xlabel("θ, °")
axes[1, 1].set_ylabel("Count")

# Plot correct detections
axes[2, 0].hist(correct_az, bins=azimuth_bins, edgecolor='black', alpha=0.7)
axes[2, 0].set_title("Correct Detections")
axes[2, 0].set_xlabel("φ,°")
axes[2, 0].set_ylabel("Count")

axes[2, 1].hist(correct_inc, bins=inclination_bins, edgecolor='black', alpha=0.7)
axes[2, 1].set_title("Correct Detections")
axes[2, 1].set_xlabel("θ, °")
axes[2, 1].set_ylabel("Count")

# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
plt.savefig("./accuracy_results.png", dpi=300)
plt.close()  

