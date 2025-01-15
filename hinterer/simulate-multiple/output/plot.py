import numpy as np
import matplotlib.pyplot as plt
from fitting_results import *

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
blob_index = np.arange(len(x_tru))

# Handle outliers

x_err, y_err, inc_err, az_err, x_est, y_est, inc_est, az_est, x_tru, y_tru, inc_tru, az_tru, blob_index = remove_outliers(x_err, y_err, inc_err, az_err, x_est, y_est, inc_est, az_est, x_tru, y_tru, inc_tru, az_tru, blob_index)
y_err, x_err, inc_err, az_err, x_est, y_est, inc_est, az_est, x_tru, y_tru, inc_tru, az_tru, blob_index = remove_outliers(y_err, x_err, inc_err, az_err, x_est, y_est, inc_est, az_est, x_tru, y_tru, inc_tru, az_tru, blob_index)

# Compute errors parallel/perpendicular to azimuth
para_err = x_err*np.sin(az_tru) + y_err*np.cos(az_tru)
perp_err = x_err*np.cos(az_tru) - y_err*np.sin(az_tru)

# Plot
fig, axs = plt.subplots(2, 2, figsize=(18, 10))

# Scatter plots
axs[0, 0].scatter(blob_index, para_err, label='hinterer')
axs[0, 0].axhline(0, color='black', linestyle='dashed', linewidth=1)
axs[0, 0].set_xlabel('Blob index')
axs[0, 0].set_ylabel('Δ$\parallel$, nm')
axs[0, 0].set_title("Localisation residuals, $\parallel$")
axs[0, 0].legend()

axs[0, 1].scatter(blob_index, perp_err, label='hinterer')
axs[0, 1].axhline(0, color='black', linestyle='dashed', linewidth=1)
axs[0, 1].set_xlabel('Blob index')
axs[0, 1].set_ylabel('Δ$\perp$, nm')
axs[0, 1].set_title("Localisation residuals, $\perp$")
axs[0, 1].legend()

#axs[0, 2].scatter(inc_tru_500, inc_err_500, label='500')
#axs[0, 2].scatter(inc_tru_1e10, inc_err_1e10, label='1e10')
#axs[0, 2].axhline(0, color='black', linestyle='dashed', linewidth=1)
#axs[0, 2].set_xlabel('θ, °')
#axs[0, 2].set_ylabel('Δθ, °')
#axs[0, 2].set_title("Inclination residual dependence on inclination")
#axs[0, 2].legend()
#
#axs[0, 3].scatter(az_tru_500, az_err_500, label='500')
#axs[0, 3].scatter(az_tru_1e10, az_err_1e10, label='1e10')
#axs[0, 3].axhline(0, color='black', linestyle='dashed', linewidth=1)
#axs[0, 3].set_xlabel('ϕ, °')
#axs[0, 3].set_ylabel('Δϕ, °')
#axs[0, 3].set_title("Azimuth residual dependence on azimuth")
#axs[0, 3].legend()

# Histogram of localisation errors

## Define equal bin sizes
#combined_min = min(np.min(x_err), np.min(y_err))
#combined_max = max(np.max(x_err), np.max(y_err))
#num_bins = 20
#bin_width = (combined_max - combined_min) / num_bins
#bins_xy = np.arange(combined_min, combined_max + bin_width, bin_width)

# Histogram of ∥ errors
axs[1, 0].hist(para_err, bins=20, alpha=0.7, label='hinterer')
axs[1, 0].axvline(np.mean(para_err), color='blue', linestyle='dashed', linewidth=1, label=f'μ = {np.mean(para_err):.2f}')
axs[1, 0].set_xlabel('Δ$\parallel$, nm')
axs[1, 0].set_ylabel('Frequency')
axs[1, 0].fill_betweenx(
    y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
    x1=np.mean(para_err) - np.std(para_err),
    x2=np.mean(para_err) + np.std(para_err),
    color='blue',
    alpha=0.1,
    label=f'σ = {np.std(para_err):.2f}'
)
axs[1, 0].set_title("$\parallel$ residual histogram")
axs[1, 0].legend()

# Histogram of ⟂ errors
axs[1, 1].hist(perp_err, bins=20, alpha=0.7, label='hinterer')
axs[1, 1].axvline(np.mean(perp_err), color='blue', linestyle='dashed', linewidth=1, label=f'μ = {np.mean(perp_err):.2f}')
axs[1, 1].set_xlabel('Δ$\perp$, nm')
axs[1, 1].set_ylabel('Frequency')
axs[1, 1].fill_betweenx(
    y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
    x1=np.mean(perp_err) - np.std(perp_err),
    x2=np.mean(perp_err) + np.std(perp_err),
    color='blue',
    alpha=0.1,
    label=f'σ = {np.std(perp_err):.2f}'
)
axs[1, 1].set_title("$\perp$ residual histogram")
axs[1, 1].legend()

## Histogram of inclination errors
#axs[1, 2].hist(inc_err_500, bins=20, alpha=0.7, label='500')
#axs[1, 2].hist(inc_err_1e10, bins=20, alpha=0.7, label='1e10')
#axs[1, 2].axvline(np.mean(inc_err_500), color='blue', linestyle='dashed', linewidth=1, label=f'μ = {np.mean(inc_err_500):.2f}')
#axs[1, 2].set_xlabel('Δθ, nm')
#axs[1, 2].set_ylabel('Frequency')
#axs[1, 2].fill_betweenx(
#    y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
#    x1=np.mean(inc_err_500) - np.std(inc_err_500),
#    x2=np.mean(inc_err_500) + np.std(inc_err_500),
#    color='blue',
#    alpha=0.1,
#    label=f'σ = {np.std(inc_err_500):.2f}'
#)
#axs[1, 2].axvline(np.mean(inc_err_1e10), color='orange', linestyle='dashed', linewidth=1, label=f'μ = {np.mean(inc_err_1e10):.2f}')
#axs[1, 2].fill_betweenx(
#    y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
#    x1=np.mean(inc_err_1e10) - np.std(inc_err_1e10),
#    x2=np.mean(inc_err_1e10) + np.std(inc_err_1e10),
#    color='orange',
#    alpha=0.1,
#    label=f'σ = {np.std(inc_err_1e10):.2f}'
#)
#axs[1, 2].set_title("θ residual histogram")
#axs[1, 2].legend()
#
## Histogram of azimuth errors
#axs[1, 3].hist(az_err_500, bins=20, alpha=0.7, label='500')
#axs[1, 3].hist(az_err_1e10, bins=20, alpha=0.7, label='1e10')
#axs[1, 3].axvline(np.mean(az_err_500), color='blue', linestyle='dashed', linewidth=1, label=f'μ = {np.mean(az_err_500):.2f}')
#axs[1, 3].set_xlabel('Δϕ, nm')
#axs[1, 3].set_ylabel('Frequency')
#axs[1, 3].fill_betweenx(
#    y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
#    x1=np.mean(az_err_500) - np.std(az_err_500),
#    x2=np.mean(az_err_500) + np.std(az_err_500),
#    color='blue',
#    alpha=0.1,
#    label=f'σ = {np.std(az_err_500):.2f}'
#)
#axs[1, 3].axvline(np.mean(az_err_1e10), color='orange', linestyle='dashed', linewidth=1, label=f'μ = {np.mean(az_err_1e10):.2f}')
#axs[1, 3].fill_betweenx(
#    y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
#    x1=np.mean(az_err_1e10) - np.std(az_err_1e10),
#    x2=np.mean(az_err_1e10) + np.std(az_err_1e10),
#    color='orange',
#    alpha=0.1,
#    label=f'σ = {np.std(az_err_1e10):.2f}'
#)
#axs[1, 3].set_title("ϕ residual histogram")
#axs[1, 3].legend()


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

