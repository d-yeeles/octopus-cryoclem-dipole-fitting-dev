import numpy as np
import matplotlib.pyplot as plt

# Function to count the number of detections in each file (TP, FP, FN)
def count_detections(file_path):
    """Counts the number of detections from a results file (TP, FP, FN)."""
    with open(file_path, 'r') as f:
        lines = f.readlines()[1:]  # Skip header if present
        return len(lines)  # Return the number of detections (rows)

# Define the threshold values and corresponding file patterns
radii = np.arange(200, 2000, 50)

precision_values_gaussian = []
recall_values_gaussian = []

# Process each radius threshold
for radius in radii:
    tp_file = f"results/TP_gaussian_radius{radius}.csv"
    fp_file = f"results/FP_gaussian_radius{radius}.csv"
    fn_file = f"results/FN_gaussian_radius{radius}.csv"
    
    tp_count = count_detections(tp_file)
    fp_count = count_detections(fp_file)
    fn_count = count_detections(fn_file)
    
    precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0
    recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0
    
    precision_values_gaussian.append(precision)
    recall_values_gaussian.append(recall)

precision_values_hinterer = []
recall_values_hinterer = []

for radius in radii:
    tp_file = f"results/TP_hinterer_radius{radius}.csv"
    fp_file = f"results/FP_hinterer_radius{radius}.csv"
    fn_file = f"results/FN_hinterer_radius{radius}.csv"
    
    tp_count = count_detections(tp_file)
    fp_count = count_detections(fp_file)
    fn_count = count_detections(fn_file)
    
    precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0
    recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0
    
    precision_values_hinterer.append(precision)
    recall_values_hinterer.append(recall)

# Plot Precision-Recall Curve
plt.figure(figsize=(8, 6))
plt.plot(recall_values_gaussian, precision_values_gaussian, marker="x", linestyle="-", label="Gaussian")
plt.plot(recall_values_hinterer, precision_values_hinterer, marker="x", linestyle="-", label="Hinterer")
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision-Recall Curve")
plt.legend(loc="upper right")
plt.grid()
plt.savefig("precision_recall.png", dpi=300)
plt.close()

