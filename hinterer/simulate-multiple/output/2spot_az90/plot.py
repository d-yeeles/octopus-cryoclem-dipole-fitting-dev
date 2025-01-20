import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
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

def count_outer_ring(xlist, ylist):
    upper_bound = 60
    x_outliers = [val for val in xlist if abs(val) > upper_bound]
    y_outliers = [val for val in ylist if abs(val) > upper_bound]
    return len(x_outliers)+len(y_outliers)

def remove_outer_ring(*lists):
    x = lists[0]
    upper_bound = 60
    outlier_indices = [i for i, val in enumerate(x) if abs(val) > upper_bound]
    return tuple([val for i, val in enumerate(lst) if i not in outlier_indices] for lst in lists)

def keep_outer_ring(*lists):
    x = lists[0]
    upper_bound = 60
    outlier_indices = [i for i, val in enumerate(x) if abs(val) < upper_bound]
    return tuple([val for i, val in enumerate(lst) if i not in outlier_indices] for lst in lists)

# Get default matplotlib colours
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
dblue = default_colors[0]
dyellow = default_colors[1]
dgreen = default_colors[2]
dred = default_colors[3]

# Create custom colormaps for hists
dBlues = LinearSegmentedColormap.from_list('dblue_to_white', [(1, 1, 1), dblue], N=100)
dYellows = LinearSegmentedColormap.from_list('dyellow_to_white', [(1, 1, 1), dyellow], N=100)

# Compute errors parallel/perpendicular to azimuth
hinterer.para_err = hinterer.x_err*np.cos(hinterer.az_tru) - hinterer.y_err*np.sin(hinterer.az_tru)
hinterer.perp_err = hinterer.x_err*np.sin(hinterer.az_tru) + hinterer.y_err*np.cos(hinterer.az_tru)
gaussian.para_err = gaussian.x_err*np.cos(gaussian.az_tru) - gaussian.y_err*np.sin(gaussian.az_tru)
gaussian.perp_err = gaussian.x_err*np.sin(gaussian.az_tru) + gaussian.y_err*np.cos(gaussian.az_tru)

# Convert inc, az to deg
convert_lists_to_degrees(hinterer, ['inc_err', 'az_err', 'inc_tru'])
convert_lists_to_degrees(gaussian, ['inc_err', 'az_err', 'inc_tru'])

# Account for wrapping of inc around 180 degrees:
gaussian.inc_err = np.minimum(np.abs(gaussian.inc_err), np.abs(180 - np.abs(gaussian.inc_err)))
hinterer.inc_err = np.minimum(np.abs(hinterer.inc_err), np.abs(180 - np.abs(hinterer.inc_err)))
# Account for wrapping of az around 360 degrees:
gaussian.az_err = np.minimum(np.abs(gaussian.az_err), np.abs(360 - np.abs(gaussian.az_err)))
hinterer.az_err = np.minimum(np.abs(hinterer.az_err), np.abs(360 - np.abs(hinterer.az_err)))

# Handle outliers - separate out the ring and the centre

N_outliers_gaussian = count_outer_ring(gaussian.para_err, gaussian.perp_err)
N_outliers_hinterer = count_outer_ring(hinterer.para_err, hinterer.perp_err)
frac_outliers_gaussian = round((N_outliers_gaussian/len(gaussian.para_err))*100, 1)
frac_outliers_hinterer = round((N_outliers_hinterer/len(hinterer.para_err))*100, 1)

gaussian.para_err_inner, gaussian.perp_err_inner, gaussian.inc_err_inner, gaussian.az_err_inner, gaussian.inc_tru_inner = remove_outer_ring(gaussian.para_err, gaussian.perp_err, gaussian.inc_err, gaussian.az_err, gaussian.inc_tru)
gaussian.perp_err_inner, gaussian.para_err_inner, gaussian.inc_err_inner, gaussian.az_err_inner, gaussian.inc_tru_inner = remove_outer_ring(gaussian.perp_err_inner, gaussian.para_err_inner, gaussian.inc_err_inner, gaussian.az_err_inner, gaussian.inc_tru_inner)
hinterer.para_err_inner, hinterer.perp_err_inner, hinterer.inc_err_inner, hinterer.az_err_inner, hinterer.inc_tru_inner = remove_outer_ring(hinterer.para_err, hinterer.perp_err, hinterer.inc_err, hinterer.az_err, hinterer.inc_tru)
hinterer.perp_err_inner, hinterer.para_err_inner, hinterer.inc_err_inner, hinterer.az_err_inner, hinterer.inc_tru_inner = remove_outer_ring(hinterer.perp_err_inner, hinterer.para_err_inner, hinterer.inc_err_inner, hinterer.az_err_inner, hinterer.inc_tru_inner)

gaussian.para_err_outer, gaussian.perp_err_outer, gaussian.inc_err_outer, gaussian.az_err_outer = keep_outer_ring(gaussian.para_err, gaussian.perp_err, gaussian.inc_err, gaussian.az_err)
gaussian.perp_err_outer, gaussian.para_err_outer, gaussian.inc_err_outer, gaussian.az_err_outer = keep_outer_ring(gaussian.perp_err_outer, gaussian.para_err_outer, gaussian.inc_err_outer, gaussian.az_err_outer)
hinterer.para_err_outer, hinterer.perp_err_outer, hinterer.inc_err_outer, hinterer.az_err_outer = keep_outer_ring(hinterer.para_err, hinterer.perp_err, hinterer.inc_err, hinterer.az_err)
hinterer.perp_err_outer, hinterer.para_err_outer, hinterer.inc_err_outer, hinterer.az_err_outer = keep_outer_ring(hinterer.perp_err_outer, hinterer.para_err_outer, hinterer.inc_err_outer, hinterer.az_err_outer)


# Plots - par/perp
fig, axs = plt.subplots(2, 3, figsize=(18, 10))

data1 = gaussian.para_err
data2 = gaussian.perp_err
mean1 = np.mean(data1)
std1 = np.std(data1)
mean2 = np.mean(data2)
std2 = np.std(data2)
axs[0, 0].scatter(data1, data2, s=2, color=dblue)
axs[0, 0].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
axs[0, 0].axhline(mean2, color=dred, linewidth=1)
axs[0, 0].add_patch(
    plt.Rectangle(
        (mean1 - std1, mean2 - std2),  # Bottom-left corner
        2 * std1, # Width
        2 * std2, # Height
        color=dred,
        fill=False,
        linestyle='dashed',
        linewidth=1,
        label=f'σ = ({std1:.2f}, {std2:.2f})'
    )
)
axs[0, 0].set_xlabel('Δ$\parallel$')
axs[0, 0].set_ylabel('Δ$\perp$')
axs[0, 0].set_title("Overview of localisation residuals")
axs[0, 0].set_title(f"Overview of localisation residuals\n{frac_outliers_gaussian}% in ring")
axs[0, 0].legend(loc='upper right')


data1 = gaussian.para_err_inner
data2 = gaussian.perp_err_inner
mean1 = np.mean(data1)
std1 = np.std(data1)
mean2 = np.mean(data2)
std2 = np.std(data2)
c = axs[0, 1].hist2d(data1, data2, bins=100, cmap=dBlues)
cbar = fig.colorbar(c[3], ax=axs[0, 1], label='Counts')
axs[0, 1].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
axs[0, 1].axhline(mean2, color=dred, linewidth=1)
axs[0, 1].add_patch(
    plt.Rectangle(
        (mean1 - std1, mean2 - std2),  # Bottom-left corner
        2 * std1, # Width
        2 * std2, # Height
        color=dred,
        fill=False,
        linestyle='dashed',
        linewidth=1,
        label=f'σ = ({std1:.2f}, {std2:.2f})'
    )
)
axs[0, 1].set_xlabel('Δ$\parallel$')
axs[0, 1].set_ylabel('Δ$\perp$')
axs[0, 1].set_title('Localisation residuals histogram - centre only')
axs[0, 1].legend(loc='upper right')

data1 = gaussian.para_err_outer
data2 = gaussian.perp_err_outer
mean1 = np.mean(data1)
std1 = np.std(data1)
mean2 = np.mean(data2)
std2 = np.std(data2)
c = axs[0, 2].hist2d(data1, data2, bins=100, cmap=dBlues)
cbar = fig.colorbar(c[3], ax=axs[0, 2], label='Counts')
cbar.set_ticks(np.arange(int(np.min(c[0])), int(np.max(c[0])) + 1, 1))  # Only show integers
axs[0, 2].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
axs[0, 2].axhline(mean2, color=dred, linewidth=1)
axs[0, 2].add_patch(
    plt.Rectangle(
        (mean1 - std1, mean2 - std2),  # Bottom-left corner
        2 * std1, # Width
        2 * std2, # Height
        color=dred,
        fill=False,
        linestyle='dashed',
        linewidth=1,
        label=f'σ = ({std1:.2f}, {std2:.2f})'
    )
)
axs[0, 2].set_xlabel('Δ$\parallel$')
axs[0, 2].set_ylabel('Δ$\perp$')
axs[0, 2].set_title('Localisation residuals histogram - ring only')
axs[0, 2].legend(loc='upper right')


data1 = hinterer.para_err
data2 = hinterer.perp_err
mean1 = np.mean(data1)
std1 = np.std(data1)
mean2 = np.mean(data2)
std2 = np.std(data2)
axs[1, 0].scatter(data1, data2, s=2, color=dyellow)
axs[1, 0].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
axs[1, 0].axhline(mean2, color=dred, linewidth=1)
axs[1, 0].add_patch(
    plt.Rectangle(
        (mean1 - std1, mean2 - std2),  # Bottom-left corner
        2 * std1, # Width
        2 * std2, # Height
        color=dred,
        fill=False,
        linestyle='dashed',
        linewidth=1,
        label=f'σ = ({std1:.2f}, {std2:.2f})'
    )
)
axs[1, 0].set_xlabel('Δ$\parallel$')
axs[1, 0].set_ylabel('Δ$\perp$')
axs[1, 0].set_title(f"Overview of localisation residuals\n{frac_outliers_hinterer}% in ring")
axs[1, 0].legend(loc='upper right')

data1 = hinterer.para_err_inner
data2 = hinterer.perp_err_inner
mean1 = np.mean(data1)
std1 = np.std(data1)
mean2 = np.mean(data2)
std2 = np.std(data2)
c = axs[1, 1].hist2d(data1, data2, bins=100, cmap=dYellows)
cbar = fig.colorbar(c[3], ax=axs[1, 1], label='Counts')
axs[1, 1].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
axs[1, 1].axhline(mean2, color=dred, linewidth=1)
axs[1, 1].add_patch(
    plt.Rectangle(
        (mean1 - std1, mean2 - std2),  # Bottom-left corner
        2 * std1, # Width
        2 * std2, # Height
        color=dred,
        fill=False,
        linestyle='dashed',
        linewidth=1,
        label=f'σ = ({std1:.2f}, {std2:.2f})'
    )
)
axs[1, 1].set_xlabel('Δ$\parallel$')
axs[1, 1].set_ylabel('Δ$\perp$')
axs[1, 1].set_title('Localisation residuals histogram - centre only')
axs[1, 1].legend(loc='upper right')

data1 = hinterer.para_err_outer
data2 = hinterer.perp_err_outer
mean1 = np.mean(data1)
std1 = np.std(data1)
mean2 = np.mean(data2)
std2 = np.std(data2)
c = axs[1, 2].hist2d(data1, data2, bins=100, cmap=dYellows)
cbar = fig.colorbar(c[3], ax=axs[1, 2], label='Counts')
cbar.set_ticks(np.arange(int(np.min(c[0])), int(np.max(c[0])) + 1, 1))  # Only show integers
axs[1, 2].axvline(mean1, color=dred, linewidth=1, label=f'μ = ({mean1:.2f}, {mean2:.2f})')
axs[1, 2].axhline(mean2, color=dred, linewidth=1)
axs[1, 2].add_patch(
    plt.Rectangle(
        (mean1 - std1, mean2 - std2),  # Bottom-left corner
        2 * std1, # Width
        2 * std2, # Height
        color=dred,
        fill=False,
        linestyle='dashed',
        linewidth=1,
        label=f'σ = ({std1:.2f}, {std2:.2f})'
    )
)
axs[1, 2].set_xlabel('Δ$\parallel$')
axs[1, 2].set_ylabel('Δ$\perp$')
axs[1, 2].set_title('Localisation residuals histogram - ring only')
axs[1, 2].legend(loc='upper right')

## Share axis scales
#x_limits = axs[0, 0].get_xlim()
#y_limits = axs[0, 0].get_ylim()
#axs[1, 0].set_xlim(x_limits)
#axs[1, 0].set_ylim(y_limits)

plt.tight_layout()








# Plots - inc/az gaussian
fig, axs = plt.subplots(2, 3, figsize=(18, 10))

# Scatter plots
axs[0, 0].scatter(gaussian.inc_tru_inner, gaussian.para_err_inner, s=10, alpha=0.1, color=dblue, label='gaussian')
axs[0, 0].axhline(0, color=dred, linewidth=2)
axs[0, 0].set_xlabel('θ, °')
axs[0, 0].set_ylabel('Δ$\parallel$, nm')
axs[0, 0].set_title("Localisation residuals, parallel - centre only")
axs[0, 0].legend(loc='upper right')

axs[0, 1].scatter(gaussian.inc_tru_inner, gaussian.perp_err_inner, s=10, alpha=0.1, color=dblue, label='gaussian')
axs[0, 1].axhline(0, color=dred, linewidth=2)
axs[0, 1].set_xlabel('θ, °')
axs[0, 1].set_ylabel('Δ$\perp$, nm')
axs[0, 1].set_title("Localisation residuals, perpendicular - centre only")
axs[0, 1].legend(loc='upper right')

axs[0, 2].scatter(gaussian.inc_tru_inner, gaussian.inc_err_inner, s=10, alpha=0.1, color=dblue, label='gaussian')
axs[0, 2].axhline(0, color=dred, linewidth=2)
axs[0, 2].set_xlabel('θ, °')
axs[0, 2].set_ylabel('Δθ, °')
axs[0, 2].set_title("Localisation residuals, polar - centre only")
axs[0, 2].legend(loc='upper right')

# Histogram of ∥ errors
data = gaussian.para_err_inner
mean = np.mean(data)
std = np.std(data)
axs[1, 0].hist(data, bins=20, alpha=0.7, color=dblue, label='gaussian')
axs[1, 0].axvline(mean, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[1, 0].fill_betweenx(
    y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dblue,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[1, 0].set_xlabel('Δ$\parallel$, nm')
axs[1, 0].set_ylabel('Frequency')
axs[1, 0].set_title("Parallel residual histogram - centre only")
axs[1, 0].legend(loc='upper right')

# Histogram of ⟂ errors
data = gaussian.perp_err_inner
mean = np.mean(data)
std = np.std(data)
axs[1, 1].hist(data, bins=20, alpha=0.7, color=dblue, label='gaussian')
axs[1, 1].axvline(mean, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[1, 1].fill_betweenx(
    y=[0, axs[1, 1].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dblue,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[1, 1].set_xlabel('Δ$\perp$, nm')
axs[1, 1].set_ylabel('Frequency')
axs[1, 1].set_title("Perpendicular residual histogram - centre only")
axs[1, 1].legend(loc='upper right')

# Histogram of inc errors
data = gaussian.inc_err_inner
mean = np.mean(data)
std = np.std(data)
axs[1, 2].hist(data, bins=20, alpha=0.7, color=dblue, label='gaussian')
axs[1, 2].axvline(mean, color=dblue, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[1, 2].fill_betweenx(
    y=[0, axs[1, 2].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dblue,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[1, 2].set_xlabel('Δθ, °')
axs[1, 2].set_ylabel('Frequency')
axs[1, 2].set_title("Polar residual histogram, centre only")
axs[1, 2].legend(loc='upper right')

plt.tight_layout()




# Plots - inc/az hinterer
fig, axs = plt.subplots(2, 3, figsize=(18, 10))

# Scatter plots
axs[0, 0].scatter(hinterer.inc_tru_inner, hinterer.para_err_inner, s=10, alpha=0.1, color=dyellow, label='hinterer')
axs[0, 0].axhline(0, color=dred, linewidth=2)
axs[0, 0].set_xlabel('θ, °')
axs[0, 0].set_ylabel('Δ$\parallel$, nm')
axs[0, 0].set_title("Localisation residuals, parallel - centre only")
axs[0, 0].legend(loc='upper right')

axs[0, 1].scatter(hinterer.inc_tru_inner, hinterer.perp_err_inner, s=10, alpha=0.1, color=dyellow, label='hinterer')
axs[0, 1].axhline(0, color=dred, linewidth=2)
axs[0, 1].set_xlabel('θ, °')
axs[0, 1].set_ylabel('Δ$\perp$, nm')
axs[0, 1].set_title("Localisation residuals, perpendicular - centre only")
axs[0, 1].legend(loc='upper right')

axs[0, 2].scatter(hinterer.inc_tru_inner, hinterer.inc_err_inner, s=10, alpha=0.1, color=dyellow, label='hinterer')
axs[0, 2].axhline(0, color=dred, linewidth=2)
axs[0, 2].set_xlabel('θ, °')
axs[0, 2].set_ylabel('Δθ, °')
axs[0, 2].set_title("Localisation residuals, polar - centre only")
axs[0, 2].legend(loc='upper right')

# Histogram of ∥ errors
data = hinterer.para_err_inner
mean = np.mean(data)
std = np.std(data)
axs[1, 0].hist(data, bins=20, alpha=0.7, color=dyellow, label='hinterer')
axs[1, 0].axvline(mean, color=dyellow, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[1, 0].fill_betweenx(
    y=[0, axs[1, 0].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dyellow,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[1, 0].set_xlabel('Δ$\parallel$, nm')
axs[1, 0].set_ylabel('Frequency')
axs[1, 0].set_title("Parallel residual histogram - centre only")
axs[1, 0].legend(loc='upper right')

# Histogram of ⟂ errors
data = hinterer.perp_err_inner
mean = np.mean(data)
std = np.std(data)
axs[1, 1].hist(data, bins=20, alpha=0.7, color=dyellow, label='hinterer')
axs[1, 1].axvline(mean, color=dyellow, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[1, 1].fill_betweenx(
    y=[0, axs[1, 1].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dyellow,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[1, 1].set_xlabel('Δ$\perp$, nm')
axs[1, 1].set_ylabel('Frequency')
axs[1, 1].set_title("Perpendicular residual histogram - centre only")
axs[1, 1].legend(loc='upper right')

# Histogram of inc errors
data = hinterer.inc_err_inner
mean = np.mean(data)
std = np.std(data)
axs[1, 2].hist(data, bins=20, alpha=0.7, color=dyellow, label='hinterer')
axs[1, 2].axvline(mean, color=dyellow, linestyle='dashed', linewidth=1, label=f'μ = {mean:.2f}')
axs[1, 2].fill_betweenx(
    y=[0, axs[1, 2].get_ylim()[1]],  # Cover full y-axis range
    x1=mean - std,
    x2=mean + std,
    color=dyellow,
    alpha=0.1,
    label=f'σ = {std:.2f}'
)
axs[1, 2].set_xlabel('Δθ, °')
axs[1, 2].set_ylabel('Frequency')
axs[1, 2].set_title("Polar residual histogram, centre only")
axs[1, 2].legend(loc='upper right')

plt.tight_layout()




plt.show()




plt.show()

