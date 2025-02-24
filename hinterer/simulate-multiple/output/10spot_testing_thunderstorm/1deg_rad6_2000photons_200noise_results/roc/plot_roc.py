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




tpr_values = []  # List to store True Positive Rates (TPR)
fpr_values = []  # List to store False Positive Rates (FPR)

# Process each radius threshold
for radius in radii:
    # Define file paths based on threshold
    tp_file = f"results/TP_gaussian_radius{radius}.csv"
    fp_file = f"results/FP_gaussian_radius{radius}.csv"
    fn_file = f"results/FN_gaussian_radius{radius}.csv"

    # Count the number of detections (TP, FP, FN)
    tp_count = count_detections(tp_file)
    fp_count = count_detections(fp_file)
    fn_count = count_detections(fn_file)

    # Compute True Positive Rate (TPR) and False Positive Rate (FPR)
    tpr = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0
    fpr = fp_count / (fp_count + tp_count) if (fp_count + tp_count) > 0 else 0  # Approximation without TN

    # Store TPR and FPR values for the current threshold
    tpr_values.append(tpr)
    fpr_values.append(fpr)

# Compute AUC (Area Under Curve) using the trapezoidal rule
auc_score = np.trapz(tpr_values, fpr_values)

fpr_values_gaussian = fpr_values
tpr_values_gaussian = tpr_values
auc_score_gaussian = auc_score




tpr_values = []  # List to store True Positive Rates (TPR)
fpr_values = []  # List to store False Positive Rates (FPR)

# Process each radius threshold
for radius in radii:
    # Define file paths based on threshold
    tp_file = f"results/TP_hinterer_radius{radius}.csv"
    fp_file = f"results/FP_hinterer_radius{radius}.csv"
    fn_file = f"results/FN_hinterer_radius{radius}.csv"

    # Count the number of detections (TP, FP, FN)
    tp_count = count_detections(tp_file)
    fp_count = count_detections(fp_file)
    fn_count = count_detections(fn_file)

    # Compute True Positive Rate (TPR) and False Positive Rate (FPR)
    tpr = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0
    fpr = fp_count / (fp_count + tp_count) if (fp_count + tp_count) > 0 else 0  # Approximation without TN

    # Store TPR and FPR values for the current threshold
    tpr_values.append(tpr)
    fpr_values.append(fpr)

# Compute AUC (Area Under Curve) using the trapezoidal rule
auc_score = np.trapz(tpr_values, fpr_values)

fpr_values_hinterer = fpr_values
tpr_values_hinterer = tpr_values
auc_score_hinterer = auc_score





# Plot ROC Curve
plt.figure(figsize=(8, 6))
plt.plot(fpr_values_gaussian, tpr_values_gaussian, marker="x", linestyle="-", label=f"Gaussian (AUC = {auc_score_gaussian:.2f})")
plt.plot(fpr_values_hinterer, tpr_values_hinterer, marker="x", linestyle="-", label=f"Hinterer (AUC = {auc_score_hinterer:.2f})")
plt.plot([0, 1], [0, 1], linestyle="--", color="gray")  # Diagonal line (random guess)
plt.xlabel("False Positive Rate (FPR)")
plt.ylabel("True Positive Rate (TPR)")
plt.title("Threshold = radius within which estimate and ground truth must be")
plt.legend(loc="lower right")
plt.grid()
#plt.show()
plt.savefig(f"roc.png", dpi=300)
plt.close()


