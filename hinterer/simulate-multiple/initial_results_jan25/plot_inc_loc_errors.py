import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from results import *

def remove_outliers(*lists):

    x = lists[0]

    q1, q3 = np.percentile(x, [25, 75])
    iqr = q3 - q1
    lower_bound = q1 - 5 * iqr
    upper_bound = q3 + 5 * iqr
    outlier_indices = [i for i, val in enumerate(x) if val < lower_bound or val > upper_bound]

    print('number of outliers excluded:', len(outlier_indices))
    
    # Remove outliers from all lists
    return tuple([val for i, val in enumerate(lst) if i not in outlier_indices] for lst in lists)

# Account for symmetry around 90 degrees:
# if it should be 0 but it estimated 80, this is only really 10 deg out because 0=90
#inc_err_wrapped_500 = np.minimum(np.abs(inc_err_500), np.abs(180 - np.abs(inc_err_500)))
#inc_err_wrapped_1e10 = np.minimum(np.abs(inc_err_1e10), np.abs(180 - np.abs(inc_err_1e10)))

# Make list of blob index
blob_index = np.arange(50)

# Handle outliers

x_err_500_optimised, y_err_500_optimised, inc_err_500_optimised, az_err_500_optimised, inc_tru_500_optimised, az_tru_500_optimised, blob_index_500_optimised = remove_outliers(x_err_500_optimised, y_err_500_optimised, inc_err_500_optimised, az_err_500_optimised, inc_tru_500_optimised, az_tru_500_optimised, blob_index)
y_err_500_optimised, x_err_500_optimised, inc_err_500_optimised, az_err_500_optimised, inc_tru_500_optimised, az_tru_500_optimised, blob_index_500_optimised = remove_outliers(y_err_500_optimised, x_err_500_optimised, inc_err_500_optimised, az_err_500_optimised, inc_tru_500_optimised, az_tru_500_optimised, blob_index_500_optimised)

x_err_500_given, y_err_500_given, inc_err_500_given, az_err_500_given, inc_tru_500_given, az_tru_500_given, blob_index_500_given = remove_outliers(x_err_500_given, y_err_500_given, inc_err_500_given, az_err_500_given, inc_tru_500_given, az_tru_500_given, blob_index)
y_err_500_given, x_err_500_given, inc_err_500_given, az_err_500_given, inc_tru_500_given, az_tru_500_given, blob_index_500_given = remove_outliers(y_err_500_given, x_err_500_given, inc_err_500_given, az_err_500_given, inc_tru_500_given, az_tru_500_given, blob_index_500_given)

x_err_500_zero, y_err_500_zero, inc_err_500_zero, az_err_500_zero, inc_tru_500_zero, az_tru_500_zero, blob_index_500_zero = remove_outliers(x_err_500_zero, y_err_500_zero, inc_err_500_zero, az_err_500_zero, inc_tru_500_zero, az_tru_500_zero, blob_index)
y_err_500_zero, x_err_500_zero, inc_err_500_zero, az_err_500_zero, inc_tru_500_zero, az_tru_500_zero, blob_index_500_zero = remove_outliers(y_err_500_zero, x_err_500_zero, inc_err_500_zero, az_err_500_zero, inc_tru_500_zero, az_tru_500_zero, blob_index_500_zero)

x_err_1e10_optimised, y_err_1e10_optimised, inc_err_1e10_optimised, az_err_1e10_optimised, inc_tru_1e10_optimised, az_tru_1e10_optimised, blob_index_1e10_optimised = remove_outliers(x_err_1e10_optimised, y_err_1e10_optimised, inc_err_1e10_optimised, az_err_1e10_optimised, inc_tru_1e10_optimised, az_tru_1e10_optimised, blob_index)
y_err_1e10_optimised, x_err_1e10_optimised, inc_err_1e10_optimised, az_err_1e10_optimised, inc_tru_1e10_optimised, az_tru_1e10_optimised, blob_index_1e10_optimised = remove_outliers(y_err_1e10_optimised, x_err_1e10_optimised, inc_err_1e10_optimised, az_err_1e10_optimised, inc_tru_1e10_optimised, az_tru_1e10_optimised, blob_index_1e10_optimised)

x_err_1e10_given, y_err_1e10_given, inc_err_1e10_given, az_err_1e10_given, inc_tru_1e10_given, az_tru_1e10_given, blob_index_1e10_given = remove_outliers(x_err_1e10_given, y_err_1e10_given, inc_err_1e10_given, az_err_1e10_given, inc_tru_1e10_given, az_tru_1e10_given, blob_index)
y_err_1e10_given, x_err_1e10_given, inc_err_1e10_given, az_err_1e10_given, inc_tru_1e10_given, az_tru_1e10_given, blob_index_1e10_given = remove_outliers(y_err_1e10_given, x_err_1e10_given, inc_err_1e10_given, az_err_1e10_given, inc_tru_1e10_given, az_tru_1e10_given, blob_index_1e10_given)

x_err_1e10_zero, y_err_1e10_zero, inc_err_1e10_zero, az_err_1e10_zero, inc_tru_1e10_zero, az_tru_1e10_zero, blob_index_1e10_zero = remove_outliers(x_err_1e10_zero, y_err_1e10_zero, inc_err_1e10_zero, az_err_1e10_zero, inc_tru_1e10_zero, az_tru_1e10_zero, blob_index)
y_err_1e10_zero, x_err_1e10_zero, inc_err_1e10_zero, az_err_1e10_zero, inc_tru_1e10_zero, az_tru_1e10_zero, blob_index_1e10_zero = remove_outliers(y_err_1e10_zero, x_err_1e10_zero, inc_err_1e10_zero, az_err_1e10_zero, inc_tru_1e10_zero, az_tru_1e10_zero, blob_index_1e10_zero)

# Mean etc.
#print('mean x:', round(np.mean(x_err), 2))
#print('std x:', round(np.std(x_err), 2))
#print('mean y:', round(np.mean(y_err), 2))
#print('std y:', round(np.std(y_err), 2))

# Plot these ones
blob_index_500 = blob_index_500_optimised
x_err_500 = x_err_500_optimised
y_err_500 = y_err_500_optimised
inc_err_500 = inc_err_500_optimised
az_err_500 = az_err_500_optimised
inc_tru_500 = inc_tru_500_optimised
az_tru_500 = az_tru_500_optimised

blob_index_1e10 = blob_index_1e10_optimised
x_err_1e10 = x_err_1e10_optimised
y_err_1e10 = y_err_1e10_optimised
inc_err_1e10 = inc_err_1e10_optimised
az_err_1e10 = az_err_1e10_optimised
inc_tru_1e10 = inc_tru_1e10_optimised
az_tru_1e10 = az_tru_1e10_optimised

# Compute errors parallel/perpendicular to azimuth
para_err_500 = x_err_500*np.sin(az_tru_500) + y_err_500*np.cos(az_tru_500)
perp_err_500 = x_err_500*np.cos(az_tru_500) - y_err_500*np.sin(az_tru_500)
para_err_1e10 = x_err_1e10*np.sin(az_tru_1e10) + y_err_1e10*np.cos(az_tru_1e10)
perp_err_1e10 = x_err_1e10*np.cos(az_tru_1e10) - y_err_1e10*np.sin(az_tru_1e10)

# Plot
fig, axs = plt.subplots(2, 3, figsize=(18, 10))

# Scatter plots
axs[0, 0].scatter(blob_index_500, para_err_500, label='500')
axs[0, 0].scatter(blob_index_1e10, para_err_1e10, label='1e10')
axs[0, 0].axhline(0, color='red', linewidth=2)
axs[0, 0].set_xlabel('Blob index')
axs[0, 0].set_ylabel('Δ||, nm')
axs[0, 0].set_title("Localisation residuals, ||")
axs[0, 0].legend(loc='upper right')

axs[0, 1].scatter(blob_index_500, perp_err_500, label='500')
axs[0, 1].scatter(blob_index_1e10, perp_err_1e10, label='1e10')
axs[0, 1].axhline(0, color='red', linewidth=2)
axs[0, 1].set_xlabel('Blob index')
axs[0, 1].set_ylabel('Δ$\perp$, nm')
axs[0, 1].set_title("Localisation residuals, $\perp$")
axs[0, 1].legend(loc='upper right')

axs[0, 2].scatter(perp_err_500, para_err_500, s=10, label='500')
axs[0, 2].scatter(perp_err_1e10, para_err_1e10, s=10, label='1e10')
axs[0, 2].scatter(0, 0, color='red', s=10)
axs[0, 2].axhline(0, color='red', linewidth=2)
axs[0, 2].axvline(0, color='red', linewidth=2)
axs[0, 2].set_xlabel('Δ$\perp$, nm')
axs[0, 2].set_ylabel('Δ||, nm')
axs[0, 2].set_title("Localisation residuals")
axs[0, 2].legend(loc='upper right')

# Histogram of localisation errors

## Define equal bin sizes
#combined_min = min(np.min(x_err), np.min(y_err))
#combined_max = max(np.max(x_err), np.max(y_err))
#num_bins = 20
#bin_width = (combined_max - combined_min) / num_bins
#bins_xy = np.arange(combined_min, combined_max + bin_width, bin_width)

# Histogram of || errors
axs[1, 0].hist(para_err_500, bins=20, alpha=0.7, label='500')
axs[1, 0].set_xlabel('Δ||, nm')
axs[1, 0].set_ylabel('Frequency')
axs[1, 0].axvline(np.mean(para_err_500), color='blue', linestyle='dashed', linewidth=1, label=f'μ = {np.mean(para_err_500):.2f}')
axs[1, 0].fill_betweenx(
    y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
    x1=np.mean(para_err_500) - np.std(para_err_500),
    x2=np.mean(para_err_500) + np.std(para_err_500),
    color='blue',
    alpha=0.1,
    label=f'σ = {np.std(para_err_500):.2f}'
)
axs[1, 0].hist(para_err_1e10, bins=20, alpha=0.7, label='1e10')
axs[1, 0].axvline(np.mean(para_err_1e10), color='orange', linestyle='dashed', linewidth=1, label=f'μ = {np.mean(para_err_1e10):.2f}')
axs[1, 0].fill_betweenx(
    y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
    x1=np.mean(para_err_1e10) - np.std(para_err_1e10),
    x2=np.mean(para_err_1e10) + np.std(para_err_1e10),
    color='orange',
    alpha=0.1,
    label=f'σ = {np.std(para_err_1e10):.2f}'
)
axs[1, 0].set_title("|| residual histogram")
axs[1, 0].legend(loc='upper right', ncol=2)

# Histogram of ⟂ errors
axs[1, 1].hist(perp_err_500, bins=20, alpha=0.7, label='500')
axs[1, 1].set_xlabel('Δ$\perp$, nm')
axs[1, 1].set_ylabel('Frequency')
axs[1, 1].axvline(np.mean(perp_err_500), color='blue', linestyle='dashed', linewidth=1, label=f'μ = {np.mean(perp_err_500):.2f}')
axs[1, 1].fill_betweenx(
    y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
    x1=np.mean(perp_err_500) - np.std(perp_err_500),
    x2=np.mean(perp_err_500) + np.std(perp_err_500),
    color='blue',
    alpha=0.1,
    label=f'σ = {np.std(perp_err_500):.2f}'
)
axs[1, 1].hist(perp_err_1e10, bins=20, alpha=0.7, label='1e10')
axs[1, 1].axvline(np.mean(perp_err_1e10), color='orange', linestyle='dashed', linewidth=1, label=f'μ = {np.mean(perp_err_1e10):.2f}')
axs[1, 1].fill_betweenx(
    y=[0, axs[1, 1].get_ylim()[1]],  # Cover full y-axis range
    x1=np.mean(perp_err_1e10) - np.std(perp_err_1e10),
    x2=np.mean(perp_err_1e10) + np.std(perp_err_1e10),
    color='orange',
    alpha=0.1,
    label=f'σ = {np.std(perp_err_1e10):.2f}'
)
axs[1, 1].set_title("$\perp$ residual histogram")
axs[1, 1].legend(loc='upper right', ncol=2)

# Parameter space coverage
axs[1, 2].scatter(inc_tru_500, az_tru_500, color='black')
axs[1, 2].set_xlabel('θ, °')
axs[1, 2].set_ylabel('ϕ, °')
axs[1, 2].set_title("Orientation parameter space coverage")



# Synchronise axes
x_min = min(axs[1, 0].get_xlim()[0], axs[1, 1].get_xlim()[0])
x_max = max(axs[1, 0].get_xlim()[1], axs[1, 1].get_xlim()[1])
y_min = min(axs[1, 0].get_ylim()[0], axs[1, 1].get_ylim()[0])
y_max = max(axs[1, 0].get_ylim()[1], axs[1, 1].get_ylim()[1])
axs[1, 0].set_xlim(x_min, x_max)
axs[1, 1].set_xlim(x_min, x_max)
axs[1, 0].set_ylim(y_min, y_max)
axs[1, 1].set_ylim(y_min, y_max)

plt.tight_layout()
plt.show()

