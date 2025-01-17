import numpy as np
import matplotlib.pyplot as plt
import fitting_results_hinterer as hinterer
import fitting_results_gaussian as gaussian

def convert_lists_to_degrees(module, variable_names):
    for name in variable_names:
        if name in vars(module):
            vars(module)[name] = [x * 180 / np.pi for x in vars(module)[name]]

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

def remove_outer_ring(*lists):
    x = lists[0]
    upper_bound = 60
    outlier_indices = [i for i, val in enumerate(x) if abs(val) > upper_bound]
    return tuple([val for i, val in enumerate(lst) if i not in outlier_indices] for lst in lists)

# Get default matplotlib colours
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
dblue = default_colors[0]
dyellow = default_colors[1]
dgreen = default_colors[2]
dred = default_colors[3]



## Make list of blob index
#blob_index = np.arange(len(x_tru))

# Handle outliers
#hinterer.x_err, hinterer.y_err, hinterer.inc_err, hinterer.az_err, hinterer.x_est, hinterer.y_est, hinterer.inc_est, hinterer.az_est, hinterer.x_tru, hinterer.y_tru, hinterer.inc_tru, hinterer.az_tru = remove_outliers(hinterer.x_err, hinterer.y_err, hinterer.inc_err, hinterer.az_err, hinterer.x_est, hinterer.y_est, hinterer.inc_est, hinterer.az_est, hinterer.x_tru, hinterer.y_tru, hinterer.inc_tru, hinterer.az_tru)
#hinterer.y_err, hinterer.x_err, hinterer.inc_err, hinterer.az_err, hinterer.x_est, hinterer.y_est, hinterer.inc_est, hinterer.az_est, hinterer.x_tru, hinterer.y_tru, hinterer.inc_tru, hinterer.az_tru = remove_outliers(hinterer.y_err, hinterer.x_err, hinterer.inc_err, hinterer.az_err, hinterer.x_est, hinterer.y_est, hinterer.inc_est, hinterer.az_est, hinterer.x_tru, hinterer.y_tru, hinterer.inc_tru, hinterer.az_tru)
#gaussian.x_err, gaussian.y_err, gaussian.inc_err, gaussian.az_err, gaussian.x_est, gaussian.y_est, gaussian.inc_est, gaussian.az_est, gaussian.x_tru, gaussian.y_tru, gaussian.inc_tru, gaussian.az_tru = remove_outliers(gaussian.x_err, gaussian.y_err, gaussian.inc_err, gaussian.az_err, gaussian.x_est, gaussian.y_est, gaussian.inc_est, gaussian.az_est, gaussian.x_tru, gaussian.y_tru, gaussian.inc_tru, gaussian.az_tru)
#gaussian.y_err, gaussian.x_err, gaussian.inc_err, gaussian.az_err, gaussian.x_est, gaussian.y_est, gaussian.inc_est, gaussian.az_est, gaussian.x_tru, gaussian.y_tru, gaussian.inc_tru, gaussian.az_tru = remove_outliers(gaussian.y_err, gaussian.x_err, gaussian.inc_err, gaussian.az_err, gaussian.x_est, gaussian.y_est, gaussian.inc_est, gaussian.az_est, gaussian.x_tru, gaussian.y_tru, gaussian.inc_tru, gaussian.az_tru)

gaussian.x_err, gaussian.y_err, gaussian.inc_err, gaussian.az_err, gaussian.x_est, gaussian.y_est, gaussian.inc_est, gaussian.az_est, gaussian.x_tru, gaussian.y_tru, gaussian.inc_tru, gaussian.az_tru = remove_outer_ring(gaussian.x_err, gaussian.y_err, gaussian.inc_err, gaussian.az_err, gaussian.x_est, gaussian.y_est, gaussian.inc_est, gaussian.az_est, gaussian.x_tru, gaussian.y_tru, gaussian.inc_tru, gaussian.az_tru)
gaussian.y_err, gaussian.x_err, gaussian.inc_err, gaussian.az_err, gaussian.x_est, gaussian.y_est, gaussian.inc_est, gaussian.az_est, gaussian.x_tru, gaussian.y_tru, gaussian.inc_tru, gaussian.az_tru = remove_outer_ring(gaussian.y_err, gaussian.x_err, gaussian.inc_err, gaussian.az_err, gaussian.x_est, gaussian.y_est, gaussian.inc_est, gaussian.az_est, gaussian.x_tru, gaussian.y_tru, gaussian.inc_tru, gaussian.az_tru)
hinterer.x_err, hinterer.y_err, hinterer.inc_err, hinterer.az_err, hinterer.x_est, hinterer.y_est, hinterer.inc_est, hinterer.az_est, hinterer.x_tru, hinterer.y_tru, hinterer.inc_tru, hinterer.az_tru = remove_outer_ring(hinterer.x_err, hinterer.y_err, hinterer.inc_err, hinterer.az_err, hinterer.x_est, hinterer.y_est, hinterer.inc_est, hinterer.az_est, hinterer.x_tru, hinterer.y_tru, hinterer.inc_tru, hinterer.az_tru)
hinterer.y_err, hinterer.x_err, hinterer.inc_err, hinterer.az_err, hinterer.x_est, hinterer.y_est, hinterer.inc_est, hinterer.az_est, hinterer.x_tru, hinterer.y_tru, hinterer.inc_tru, hinterer.az_tru = remove_outer_ring(hinterer.y_err, hinterer.x_err, hinterer.inc_err, hinterer.az_err, hinterer.x_est, hinterer.y_est, hinterer.inc_est, hinterer.az_est, hinterer.x_tru, hinterer.y_tru, hinterer.inc_tru, hinterer.az_tru)

# Compute errors parallel/perpendicular to azimuth
hinterer.para_err = hinterer.x_err*np.cos(hinterer.az_tru) - hinterer.y_err*np.sin(hinterer.az_tru)
hinterer.perp_err = hinterer.x_err*np.sin(hinterer.az_tru) + hinterer.y_err*np.cos(hinterer.az_tru)
gaussian.para_err = gaussian.x_err*np.cos(gaussian.az_tru) - gaussian.y_err*np.sin(gaussian.az_tru)
gaussian.perp_err = gaussian.x_err*np.sin(gaussian.az_tru) + gaussian.y_err*np.cos(gaussian.az_tru)

# Convert inc, az to deg
convert_lists_to_degrees(hinterer, ['inc_tru', 'inc_est', 'inc_err', 'az_tru', 'az_est', 'az_err'])
convert_lists_to_degrees(gaussian, ['inc_tru', 'inc_est', 'inc_err', 'az_tru', 'az_est', 'az_err'])

# Account for wrapping of inc around 180 degrees:
gaussian.inc_err = np.minimum(np.abs(gaussian.inc_err), np.abs(180 - np.abs(gaussian.inc_err)))
hinterer.inc_err = np.minimum(np.abs(hinterer.inc_err), np.abs(180 - np.abs(hinterer.inc_err)))
# Account for wrapping of az around 360 degrees:
gaussian.az_err = np.minimum(np.abs(gaussian.az_err), np.abs(360 - np.abs(gaussian.az_err)))
hinterer.az_err = np.minimum(np.abs(hinterer.az_err), np.abs(360 - np.abs(hinterer.az_err)))


# Initial plot
fig, axs = plt.subplots(1, 2, figsize=(9, 5))

axs[0].scatter(gaussian.perp_err, gaussian.para_err, s=2, color=dblue)
axs[0].scatter(0, 0, color=dred, s=10)
axs[0].axhline(0, color=dred, linewidth=2)
axs[0].axvline(0, color=dred, linewidth=2)
axs[0].set_xlabel('Δ$\perp$, nm')
axs[0].set_ylabel('Δ||, nm')
axs[0].set_title("Gaussian localisation residuals")

axs[1].scatter(hinterer.perp_err, hinterer.para_err, s=2, color=dyellow)
axs[1].scatter(0, 0, color=dred, s=10)
axs[1].axhline(0, color=dred, linewidth=2)
axs[1].axvline(0, color=dred, linewidth=2)
axs[1].set_xlabel('Δ$\perp$, nm')
axs[1].set_ylabel('Δ||, nm')
axs[1].set_title("Hinterer localisation residuals")

# Share axis scales
x_limits = axs[0].get_xlim()
y_limits = axs[0].get_ylim()
axs[1].set_xlim(x_limits)
axs[1].set_ylim(y_limits)

plt.tight_layout()
#plt.show()




# # Plot - overlay
# fig, axs = plt.subplots(3, 4, figsize=(18, 10))
# 
# # Scatter plots
# axs[0, 0].scatter(gaussian.inc_tru, gaussian.para_err, s=10, alpha=0.1, label='gaussian')
# axs[0, 0].scatter(hinterer.inc_tru, hinterer.para_err, s=10, alpha=0.1, label='hinterer')
# axs[0, 0].axhline(0, color=dred, linewidth=2)
# axs[0, 0].set_xlabel('θ, °')
# axs[0, 0].set_ylabel('Δ$\parallel$, nm')
# axs[0, 0].set_title("Localisation residuals, $\parallel$")
# axs[0, 0].legend(loc='upper right')
# 
# axs[0, 1].scatter(gaussian.inc_tru, gaussian.perp_err, s=10, alpha=0.1, label='gaussian')
# axs[0, 1].scatter(hinterer.inc_tru, hinterer.perp_err, s=10, alpha=0.1, label='hinterer')
# axs[0, 1].axhline(0, color=dred, linewidth=2)
# axs[0, 1].set_xlabel('θ, °')
# axs[0, 1].set_ylabel('Δ$\perp$, nm')
# axs[0, 1].set_title("Localisation residuals, $\perp$")
# axs[0, 1].legend(loc='upper right')
# 
# axs[0, 2].scatter(gaussian.inc_tru, gaussian.inc_err, s=10, alpha=0.1, label='gaussian')
# axs[0, 2].scatter(hinterer.inc_tru, hinterer.inc_err, s=10, alpha=0.1, label='hinterer')
# axs[0, 2].axhline(0, color=dred, linewidth=2)
# axs[0, 2].set_xlabel('θ, °')
# axs[0, 2].set_ylabel('Δθ, nm')
# axs[0, 2].set_title("Localisation residuals, θ")
# axs[0, 2].legend(loc='upper right')
# 
# axs[0, 3].scatter(gaussian.inc_tru, gaussian.az_err, s=10, alpha=0.1, label='gaussian')
# axs[0, 3].scatter(hinterer.inc_tru, hinterer.az_err, s=10, alpha=0.1, label='hinterer')
# axs[0, 3].axhline(0, color=dred, linewidth=2)
# axs[0, 3].set_xlabel('θ, °')
# axs[0, 3].set_ylabel('Δφ, nm')
# axs[0, 3].set_title("Localisation residuals, φ")
# axs[0, 3].legend(loc='upper right')
# 
# 
# axs[1, 0].scatter(gaussian.az_tru, gaussian.para_err, s=10, alpha=0.1, label='gaussian')
# axs[1, 0].scatter(hinterer.az_tru, hinterer.para_err, s=10, alpha=0.1, label='hinterer')
# axs[1, 0].axhline(0, color=dred, linewidth=2)
# axs[1, 0].set_xlabel('Azimuth, °')
# axs[1, 0].set_ylabel('Δ$\parallel$, nm')
# axs[1, 0].set_title("Localisation residuals, $\parallel$")
# axs[1, 0].legend(loc='upper right')
# 
# axs[1, 1].scatter(gaussian.az_tru, gaussian.perp_err, s=10, alpha=0.1, label='gaussian')
# axs[1, 1].scatter(hinterer.az_tru, hinterer.perp_err, s=10, alpha=0.1, label='hinterer')
# axs[1, 1].axhline(0, color=dred, linewidth=2)
# axs[1, 1].set_xlabel('Azimuth, °')
# axs[1, 1].set_ylabel('Δ$\perp$, nm')
# axs[1, 1].set_title("Localisation residuals, $\perp$")
# axs[1, 1].legend(loc='upper right')
# 
# axs[1, 2].scatter(gaussian.az_tru, gaussian.inc_err, s=10, alpha=0.1, label='gaussian')
# axs[1, 2].scatter(hinterer.az_tru, hinterer.inc_err, s=10, alpha=0.1, label='hinterer')
# axs[1, 2].axhline(0, color=dred, linewidth=2)
# axs[1, 2].set_xlabel('φ, °')
# axs[1, 2].set_ylabel('Δθ, nm')
# axs[1, 2].set_title("Localisation residuals, θ")
# axs[1, 2].legend(loc='upper right')
# 
# axs[1, 3].scatter(gaussian.az_tru, gaussian.az_err, s=10, alpha=0.1, label='gaussian')
# axs[1, 3].scatter(hinterer.az_tru, hinterer.az_err, s=10, alpha=0.1, label='hinterer')
# axs[1, 3].axhline(0, color=dred, linewidth=2)
# axs[1, 3].set_xlabel('φ, °')
# axs[1, 3].set_ylabel('Δφ, nm')
# axs[1, 3].set_title("Localisation residuals, φ")
# axs[1, 3].legend(loc='upper right')
# 
# 
# # Histogram of localisation errors
# 
# ## Define equal bin sizes
# #combined_min = min(np.min(x_err), np.min(y_err))
# #combined_max = max(np.max(x_err), np.max(y_err))
# #num_bins = 20
# #bin_width = (combined_max - combined_min) / num_bins
# #bins_xy = np.arange(combined_min, combined_max + bin_width, bin_width)
# 
# # Histogram of ∥ errors
# data = gaussian.para_err
# mean = np.mean(data)
# std = np.std(data)
# axs[2, 0].hist(data, bins=20, alpha=0.7, label='gaussian')
# axs[2, 0].axvline(mean, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
# axs[2, 0].fill_betweenx(
#     y=[0, axs[2, 0].get_ylim()[1]],  # Cover full y-axis range
#     x1=mean - std,
#     x2=mean + std,
#     color=dblue,
#     alpha=0.1,
#     label=f'σ = {std:.2f}'
# )
# data = hinterer.para_err
# mean = np.mean(data)
# std = np.std(data)
# axs[2, 0].hist(data, bins=20, alpha=0.7, label='hinterer')
# axs[2, 0].axvline(mean, color=dyellow, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
# axs[2, 0].fill_betweenx(
#     y=[0, axs[2, 0].get_ylim()[1]],  # Cover full y-axis range
#     x1=mean - std,
#     x2=mean + std,
#     color=dyellow,
#     alpha=0.1,
#     label=f'σ = {std:.2f}'
# )
# axs[2, 0].set_xlabel('Δ$\parallel$, nm')
# axs[2, 0].set_ylabel('Frequency')
# axs[2, 0].set_title("$\parallel$ residual histogram")
# axs[2, 0].legend(loc='upper right', ncols=2)
# 
# # Histogram of ⟂ errors
# data = gaussian.perp_err
# mean = np.mean(data)
# std = np.std(data)
# axs[2, 1].hist(data, bins=20, alpha=0.7, label='gaussian')
# axs[2, 1].axvline(mean, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
# axs[2, 1].fill_betweenx(
#     y=[0, axs[2, 1].get_ylim()[1]],  # Cover full y-axis range
#     x1=mean - std,
#     x2=mean + std,
#     color=dblue,
#     alpha=0.1,
#     label=f'σ = {std:.2f}'
# )
# data = hinterer.para_err
# mean = np.mean(data)
# std = np.std(data)
# axs[2, 1].hist(data, bins=20, alpha=0.7, label='hinterer')
# axs[2, 1].axvline(mean, color=dyellow, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
# axs[2, 1].fill_betweenx(
#     y=[0, axs[2, 1].get_ylim()[1]],  # Cover full y-axis range
#     x1=mean - std,
#     x2=mean + std,
#     color=dyellow,
#     alpha=0.1,
#     label=f'σ = {std:.2f}'
# )
# axs[2, 1].set_xlabel('Δ$\perp$, nm')
# axs[2, 1].set_ylabel('Frequency')
# axs[2, 1].set_title("$\perp$ residual histogram")
# axs[2, 1].legend(loc='upper right', ncols=2)
# 
# # Histogram of inc errors
# data = gaussian.inc_err
# mean = np.mean(data)
# std = np.std(data)
# axs[2, 2].hist(data, bins=20, alpha=0.7, label='gaussian')
# axs[2, 2].axvline(mean, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
# axs[2, 2].fill_betweenx(
#     y=[0, axs[2, 2].get_ylim()[1]],  # Cover full y-axis range
#     x1=mean - std,
#     x2=mean + std,
#     color=dblue,
#     alpha=0.1,
#     label=f'σ = {std:.2f}'
# )
# data = hinterer.inc_err
# mean = np.mean(data)
# std = np.std(data)
# axs[2, 2].hist(data, bins=20, alpha=0.7, label='hinterer')
# axs[2, 2].axvline(mean, color=dyellow, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
# axs[2, 2].fill_betweenx(
#     y=[0, axs[2, 2].get_ylim()[1]],  # Cover full y-axis range
#     x1=mean - std,
#     x2=mean + std,
#     color=dyellow,
#     alpha=0.1,
#     label=f'σ = {std:.2f}'
# )
# axs[2, 2].set_xlabel('Δθ, nm')
# axs[2, 2].set_ylabel('Frequency')
# axs[2, 2].set_title("θ residual histogram")
# axs[2, 2].legend(loc='upper right', ncols=2)
# 
# # Histogram of az errors
# data = gaussian.az_err
# mean = np.mean(data)
# std = np.std(data)
# axs[2, 3].hist(data, bins=20, alpha=0.7, label='gaussian')
# axs[2, 3].axvline(mean, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
# axs[2, 3].fill_betweenx(
#     y=[0, axs[2, 3].get_ylim()[1]],  # Cover full y-axis range
#     x1=mean - std,
#     x2=mean + std,
#     color=dblue,
#     alpha=0.1,
#     label=f'σ = {std:.2f}'
# )
# data = hinterer.az_err
# mean = np.mean(data)
# std = np.std(data)
# axs[2, 3].hist(data, bins=20, alpha=0.7, label='hinterer')
# axs[2, 3].axvline(mean, color=dyellow, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
# axs[2, 3].fill_betweenx(
#     y=[0, axs[2, 3].get_ylim()[1]],  # Cover full y-axis range
#     x1=mean - std,
#     x2=mean + std,
#     color=dyellow,
#     alpha=0.1,
#     label=f'σ = {std:.2f}'
# )
# axs[2, 3].set_xlabel('Δφ, nm')
# axs[2, 3].set_ylabel('Frequency')
# axs[2, 3].set_title("φ residual histogram")
# axs[2, 3].legend(loc='upper right', ncols=2)
# 
# 
# # Synchronise axes
# #x_min = min(axs[1, 0].get_xlim()[0], axs[1, 1].get_xlim()[0])
# #x_max = max(axs[1, 0].get_xlim()[1], axs[1, 1].get_xlim()[1])
# #y_min = min(axs[1, 0].get_ylim()[0], axs[1, 1].get_ylim()[0])
# #y_max = max(axs[1, 0].get_ylim()[1], axs[1, 1].get_ylim()[1])
# #axs[1, 0].set_xlim(x_min, x_max)
# #axs[1, 1].set_xlim(x_min, x_max)
# #axs[1, 0].set_ylim(y_min, y_max)
# #axs[1, 1].set_ylim(y_min, y_max)
# 
# plt.tight_layout()
# plt.show()



# Plot - gaussian
fig, axs = plt.subplots(3, 4, figsize=(18, 10))

# Scatter plots
axs[0, 0].scatter(gaussian.inc_tru, gaussian.para_err, s=10, alpha=0.1, label='gaussian')
axs[0, 0].axhline(0, color=dred, linewidth=2)
axs[0, 0].set_xlabel('θ, °')
axs[0, 0].set_ylabel('Δ$\parallel$, nm')
axs[0, 0].set_title("Localisation residuals, parallel")
axs[0, 0].legend(loc='upper right')

axs[0, 1].scatter(gaussian.inc_tru, gaussian.perp_err, s=10, alpha=0.1, label='gaussian')
axs[0, 1].axhline(0, color=dred, linewidth=2)
axs[0, 1].set_xlabel('θ, °')
axs[0, 1].set_ylabel('Δ$\perp$, nm')
axs[0, 1].set_title("Localisation residuals, perpendicular")
axs[0, 1].legend(loc='upper right')

axs[0, 2].scatter(gaussian.inc_tru, gaussian.inc_err, s=10, alpha=0.1, label='gaussian')
axs[0, 2].axhline(0, color=dred, linewidth=2)
axs[0, 2].set_xlabel('θ, °')
axs[0, 2].set_ylabel('Δθ, °')
axs[0, 2].set_title("Localisation residuals, polar")
axs[0, 2].legend(loc='upper right')

axs[0, 3].axis('off')


axs[1, 0].scatter(gaussian.az_tru, gaussian.para_err, s=10, alpha=0.1, label='gaussian')
axs[1, 0].axhline(0, color=dred, linewidth=2)
axs[1, 0].set_xlabel('φ, °')
axs[1, 0].set_ylabel('Δ$\parallel$, nm')
axs[1, 0].set_title("Localisation residuals, parallel")
axs[1, 0].legend(loc='upper right')

axs[1, 1].scatter(gaussian.az_tru, gaussian.perp_err, s=10, alpha=0.1, label='gaussian')
axs[1, 1].axhline(0, color=dred, linewidth=2)
axs[1, 1].set_xlabel('φ, °')
axs[1, 1].set_ylabel('Δ$\perp$, nm')
axs[1, 1].set_title("Localisation residuals, perpendicular")
axs[1, 1].legend(loc='upper right')

axs[1, 2].axis('off')

axs[1, 3].scatter(gaussian.az_tru, gaussian.az_err, s=10, alpha=0.1, label='gaussian')
axs[1, 3].axhline(0, color=dred, linewidth=2)
axs[1, 3].set_xlabel('φ, °')
axs[1, 3].set_ylabel('Δφ, °')
axs[1, 3].set_title("Localisation residuals, azimuth")
axs[1, 3].legend(loc='upper right')


# Histogram of localisation errors

## Define equal bin sizes
#combined_min = min(np.min(x_err), np.min(y_err))
#combined_max = max(np.max(x_err), np.max(y_err))
#num_bins = 20
#bin_width = (combined_max - combined_min) / num_bins
#bins_xy = np.arange(combined_min, combined_max + bin_width, bin_width)

# Histogram of ∥ errors
data = gaussian.para_err
mean = np.mean(data)
std = np.std(data)
axs[2, 0].hist(data, bins=20, alpha=0.7, label='gaussian')
axs[2, 0].axvline(mean, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[2, 0].fill_betweenx(
    y=[0, axs[2, 0].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dblue,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[2, 0].set_xlabel('Δ$\parallel$, nm')
axs[2, 0].set_ylabel('Frequency')
axs[2, 0].set_title("Parallel residual histogram")
axs[2, 0].legend(loc='upper right')

# Histogram of ⟂ errors
data = gaussian.perp_err
mean = np.mean(data)
std = np.std(data)
axs[2, 1].hist(data, bins=20, alpha=0.7, label='gaussian')
axs[2, 1].axvline(mean, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[2, 1].fill_betweenx(
    y=[0, axs[2, 1].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dblue,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[2, 1].set_xlabel('Δ$\perp$, nm')
axs[2, 1].set_ylabel('Frequency')
axs[2, 1].set_title("Perpendicular residual histogram")
axs[2, 1].legend(loc='upper right')

# Histogram of inc errors
data = gaussian.inc_err
mean = np.mean(data)
std = np.std(data)
axs[2, 2].hist(data, bins=20, alpha=0.7, label='gaussian')
axs[2, 2].axvline(mean, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[2, 2].fill_betweenx(
    y=[0, axs[2, 2].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dblue,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[2, 2].set_xlabel('Δθ, °')
axs[2, 2].set_ylabel('Frequency')
axs[2, 2].set_title("Polar residual histogram")
axs[2, 2].legend(loc='upper right')

# Histogram of az errors
data = gaussian.az_err
mean = np.mean(data)
std = np.std(data)
axs[2, 3].hist(data, bins=20, alpha=0.7, label='gaussian')
axs[2, 3].axvline(mean, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[2, 3].fill_betweenx(
    y=[0, axs[2, 3].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dblue,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[2, 3].set_xlabel('Δφ, °')
axs[2, 3].set_ylabel('Frequency')
axs[2, 3].set_title("Azimuth residual histogram")
axs[2, 3].legend(loc='upper right')


# Synchronise axes
#x_min = min(axs[1, 0].get_xlim()[0], axs[1, 1].get_xlim()[0])
#x_max = max(axs[1, 0].get_xlim()[1], axs[1, 1].get_xlim()[1])
#y_min = min(axs[1, 0].get_ylim()[0], axs[1, 1].get_ylim()[0])
#y_max = max(axs[1, 0].get_ylim()[1], axs[1, 1].get_ylim()[1])
#axs[1, 0].set_xlim(x_min, x_max)
#axs[1, 1].set_xlim(x_min, x_max)
#axs[1, 0].set_ylim(y_min, y_max)
#axs[1, 1].set_ylim(y_min, y_max)

plt.tight_layout()



# Plot - hinterer
fig, axs = plt.subplots(3, 4, figsize=(18, 10))

# Scatter plots
axs[0, 0].scatter(hinterer.inc_tru, hinterer.para_err, s=10, alpha=0.1, color=dyellow, label='hinterer')
axs[0, 0].axhline(0, color=dred, linewidth=2)
axs[0, 0].set_xlabel('θ, °')
axs[0, 0].set_ylabel('Δ$\parallel$, nm')
axs[0, 0].set_title("Localisation residuals, parallel")
axs[0, 0].legend(loc='upper right')

axs[0, 1].scatter(hinterer.inc_tru, hinterer.perp_err, s=10, alpha=0.1, color=dyellow, label='hinterer')
axs[0, 1].axhline(0, color=dred, linewidth=2)
axs[0, 1].set_xlabel('θ, °')
axs[0, 1].set_ylabel('Δ$\perp$, nm')
axs[0, 1].set_title("Localisation residuals, perpendicular")
axs[0, 1].legend(loc='upper right')

axs[0, 2].scatter(hinterer.inc_tru, hinterer.inc_err, s=10, alpha=0.1, color=dyellow, label='hinterer')
axs[0, 2].axhline(0, color=dred, linewidth=2)
axs[0, 2].set_xlabel('θ, °')
axs[0, 2].set_ylabel('Δθ, °')
axs[0, 2].set_title("Localisation residuals, polar")
axs[0, 2].legend(loc='upper right')

axs[0, 3].axis('off')


axs[1, 0].scatter(hinterer.az_tru, hinterer.para_err, s=10, alpha=0.1, color=dyellow, label='hinterer')
axs[1, 0].axhline(0, color=dred, linewidth=2)
axs[1, 0].set_xlabel('φ, °')
axs[1, 0].set_ylabel('Δ$\parallel$, nm')
axs[1, 0].set_title("Localisation residuals, parallel")
axs[1, 0].legend(loc='upper right')

axs[1, 1].scatter(hinterer.az_tru, hinterer.perp_err, s=10, alpha=0.1, color=dyellow, label='hinterer')
axs[1, 1].axhline(0, color=dred, linewidth=2)
axs[1, 1].set_xlabel('φ, °')
axs[1, 1].set_ylabel('Δ$\perp$, nm')
axs[1, 1].set_title("Localisation residuals, perpendicular")
axs[1, 1].legend(loc='upper right')

axs[1, 2].axis('off')

axs[1, 3].scatter(hinterer.az_tru, hinterer.az_err, s=10, alpha=0.1, color=dyellow, label='hinterer')
axs[1, 3].axhline(0, color=dred, linewidth=2)
axs[1, 3].set_xlabel('φ, °')
axs[1, 3].set_ylabel('Δφ, °')
axs[1, 3].set_title("Localisation residuals, azimuth")
axs[1, 3].legend(loc='upper right')


# Histogram of localisation errors

## Define equal bin sizes
#combined_min = min(np.min(x_err), np.min(y_err))
#combined_max = max(np.max(x_err), np.max(y_err))
#num_bins = 20
#bin_width = (combined_max - combined_min) / num_bins
#bins_xy = np.arange(combined_min, combined_max + bin_width, bin_width)

# Histogram of ∥ errors
data = hinterer.para_err
mean = np.mean(data)
std = np.std(data)
axs[2, 0].hist(data, bins=20, alpha=0.7, color=dyellow, label='hinterer')
axs[2, 0].axvline(mean, color=dyellow, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[2, 0].fill_betweenx(
    y=[0, axs[2, 0].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dyellow,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[2, 0].set_xlabel('Δ$\parallel$, nm')
axs[2, 0].set_ylabel('Frequency')
axs[2, 0].set_title("Parallel residual histogram")
axs[2, 0].legend(loc='upper right')

# Histogram of ⟂ errors
data = hinterer.perp_err
mean = np.mean(data)
std = np.std(data)
axs[2, 1].hist(data, bins=20, alpha=0.7, color=dyellow, label='hinterer')
axs[2, 1].axvline(mean, color=dyellow, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[2, 1].fill_betweenx(
    y=[0, axs[2, 1].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dyellow,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[2, 1].set_xlabel('Δ$\perp$, nm')
axs[2, 1].set_ylabel('Frequency')
axs[2, 1].set_title("Perpendicular residual histogram")
axs[2, 1].legend(loc='upper right')

# Histogram of inc errors
data = hinterer.inc_err
mean = np.mean(data)
std = np.std(data)
axs[2, 2].hist(data, bins=20, alpha=0.7, color=dyellow, label='hinterer')
axs[2, 2].axvline(mean, color=dyellow, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[2, 2].fill_betweenx(
    y=[0, axs[2, 2].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dyellow,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[2, 2].set_xlabel('Δθ, °')
axs[2, 2].set_ylabel('Frequency')
axs[2, 2].set_title("Polar residual histogram")
axs[2, 2].legend(loc='upper right')

# Histogram of az errors
data = hinterer.az_err
mean = np.mean(data)
std = np.std(data)
axs[2, 3].hist(data, bins=20, alpha=0.7, color=dyellow, label='hinterer')
axs[2, 3].axvline(mean, color=dyellow, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[2, 3].fill_betweenx(
    y=[0, axs[2, 3].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dyellow,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[2, 3].set_xlabel('Δφ, °')
axs[2, 3].set_ylabel('Frequency')
axs[2, 3].set_title("Azimuth residual histogram")
axs[2, 3].legend(loc='upper right')


# Synchronise axes
#x_min = min(axs[1, 0].get_xlim()[0], axs[1, 1].get_xlim()[0])
#x_max = max(axs[1, 0].get_xlim()[1], axs[1, 1].get_xlim()[1])
#y_min = min(axs[1, 0].get_ylim()[0], axs[1, 1].get_ylim()[0])
#y_max = max(axs[1, 0].get_ylim()[1], axs[1, 1].get_ylim()[1])
#axs[1, 0].set_xlim(x_min, x_max)
#axs[1, 1].set_xlim(x_min, x_max)
#axs[1, 0].set_ylim(y_min, y_max)
#axs[1, 1].set_ylim(y_min, y_max)

plt.tight_layout()

plt.show()
