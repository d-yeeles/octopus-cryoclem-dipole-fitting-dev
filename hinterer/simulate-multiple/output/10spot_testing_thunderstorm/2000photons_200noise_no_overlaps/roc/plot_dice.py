import numpy as np
import matplotlib.pyplot as plt

# Function to count the number of detections in each file (TP, FP, FN)
def count_detections(file_path):
    """Counts the number of detections from a results file (TP, FP, FN)."""
    with open(file_path, 'r') as f:
        lines = f.readlines()[1:]  # Skip header if present
        return len(lines)  # Return the number of detections (rows)

dice_values_gaussian = []

# Process each radius threshold
for radius in radii:
    tp_file = f"results/TP_gaussian_radius{radius}.csv"
    fp_file = f"results/FP_gaussian_radius{radius}.csv"
    fn_file = f"results/FN_gaussian_radius{radius}.csv"
    
    tp_count = count_detections(tp_file)
    fp_count = count_detections(fp_file)
    fn_count = count_detections(fn_file)
    
    dice_coeff = 2*tp_count / (2*tp_count + fp_count + fn_count) if (2*tp_count + fp_count + fn_count) > 0 else 0
    
    dice_values_gaussian.append(dice_coeff)

dice_values_hinterer = []

for radius in radii:
    tp_file = f"results/TP_hinterer_radius{radius}.csv"
    fp_file = f"results/FP_hinterer_radius{radius}.csv"
    fn_file = f"results/FN_hinterer_radius{radius}.csv"
    
    tp_count = count_detections(tp_file)
    fp_count = count_detections(fp_file)
    fn_count = count_detections(fn_file)
    
    dice_coeff = 2*tp_count / (2*tp_count + fp_count + fn_count) if (2*tp_count + fp_count + fn_count) > 0 else 0
    
    dice_values_hinterer.append(dice_coeff)

# Plot
plt.figure(figsize=(8, 6))
plt.plot(recall_values_gaussian, dice_values_gaussian, marker="x", linestyle="-", label="Gaussian")
plt.plot(recall_values_hinterer, dice_values_hinterer, marker="x", linestyle="-", label="Hinterer")
plt.xlabel("Inclusion Radius, nm")
plt.ylabel("Dice, 2TP/(2TP+FP+FN)")
plt.title("Precision-Recall Curve")
plt.legend(loc="upper right")
plt.grid()
plt.savefig("dice.png", dpi=300)
plt.close()

