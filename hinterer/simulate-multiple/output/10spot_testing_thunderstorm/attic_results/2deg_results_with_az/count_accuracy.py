import pandas as pd
import numpy as np

# Load the CSV files
thunderstorm_df = pd.read_csv('thunderstorm_results.csv')
ground_truth_df = pd.read_csv('ground_truth.csv')

# Define the radius for matching points (500 nm)
radius = 500

def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

correct = []
false_positives = []
false_negatives = []

for frame in ground_truth_df['frame'].unique():

    gt_points = ground_truth_df[ground_truth_df['frame'] == frame]
    ts_points = thunderstorm_df[thunderstorm_df['frame'] == frame]
    
    # First, cycle through each TS point and compare it to every GT point
    # If there is a GT point within radius, add that GT point to 'correct'
    # If there is no GT point within radius, add that GT point to 'false negatives'

    for _, ts_row in ts_points.iterrows():
        distances = gt_points.apply(lambda row: calculate_distance(ts_row['x [nm]'], ts_row['y [nm]'], row['x_tru'], row['y_tru']), axis=1)

        if (distances <= radius).sum() == 1:
            # Save the GT row, not TS row
            correct.append(gt_points.loc[distances.idxmin()])  # Save the GT row corresponding to the nearest TS point
        elif (distances <= radius).sum() == 0:
            # Save the GT row, not TS row
            false_negatives.append(gt_points.loc[distances.idxmin()])  # Save the GT row corresponding to the nearest TS point

    # Next, cycle through each GT point and compare it to every TS point
    # If there is more than 1 TS point within radius, add that GT point to 'false positives'

    for _, gt_row in gt_points.iterrows():
        distances = ts_points.apply(lambda row: calculate_distance(gt_row['x_tru'], gt_row['y_tru'], row['x [nm]'], row['y [nm]']), axis=1)

        if (distances <= radius).sum() > 1:
            # Save the GT row, not TS row
            false_positives.append(gt_row)

pd.DataFrame(correct).to_csv('correct.csv', index=False)
pd.DataFrame(false_positives).to_csv('false_positives.csv', index=False)
pd.DataFrame(false_negatives).to_csv('false_negatives.csv', index=False)

print("Processing complete.")




# ------------------------------------------------


import matplotlib.pyplot as plt

def safe_read_csv(file):
    try:
        df = pd.read_csv(file)
        if df.empty:
            print(f"Warning: {file} is empty.")
        return df
    except pd.errors.EmptyDataError:
        print(f"Warning: {file} is empty or missing.")
        return pd.DataFrame()  # Return an empty DataFrame to prevent errors

# Read the data
ground_truth_df = safe_read_csv('ground_truth.csv')
correct_df = safe_read_csv('correct.csv')
false_positives_df = safe_read_csv('false_positives.csv')
false_negatives_df = safe_read_csv('false_negatives.csv')

# Convert radians to degrees (only if the columns exist)
for df in [ground_truth_df, correct_df, false_positives_df, false_negatives_df]:
    if 'inc_tru' in df.columns and not df.empty:
        df['inc_tru'] *= 180 / np.pi
    if 'az_tru' in df.columns and not df.empty:
        df['az_tru'] *= 180 / np.pi

# Define the azimuth bins (0 to 359, in steps of 10)
inc_bins = ground_truth_df['inc_tru'].unique()
az_bins = ground_truth_df['az_tru'].unique()

# Count occurrences of azimuth values in each bin for each category
correct_counts_inc, _ = np.histogram(correct_df['inc_tru'], bins=inc_bins, range=(0, 90))
false_positives_counts_inc, _ = np.histogram(false_positives_df['inc_tru'], bins=inc_bins, range=(0, 90))
false_negatives_counts_inc, _ = np.histogram(false_negatives_df['inc_tru'], bins=inc_bins, range=(0, 90))

# Total counts per bin (sum of all categories for that bin)
total_counts_inc = correct_counts_inc + false_positives_counts_inc + false_negatives_counts_inc

# Normalize the counts so that the sum of each bin is 1
correct_proportions_inc = (correct_counts_inc / total_counts_inc)*100
false_positives_proportions_inc = (false_positives_counts_inc / total_counts_inc)*100
false_negatives_proportions_inc = (false_negatives_counts_inc / total_counts_inc)*100




# Count occurrences of azimuth values in each bin for each category
correct_counts_az, _ = np.histogram(correct_df['az_tru'], bins=az_bins, range=(0, 360))
false_positives_counts_az, _ = np.histogram(false_positives_df['az_tru'], bins=az_bins, range=(0, 360))
false_negatives_counts_az, _ = np.histogram(false_negatives_df['az_tru'], bins=az_bins, range=(0, 360))

# Total counts per bin (sum of all categories for that bin)
total_counts_az = correct_counts_az + false_positives_counts_az + false_negatives_counts_az

# Normalize the counts so that the sum of each bin is 1
correct_proportions_az = correct_counts_az / total_counts_az
false_positives_proportions_az = false_positives_counts_az / total_counts_az
false_negatives_proportions_az = false_negatives_counts_az / total_counts_az







# Plot the results
bar_width_inc = ground_truth_df['inc_tru'].unique()[1] - ground_truth_df['inc_tru'].unique()[0]
bar_width_az = ground_truth_df['az_tru'].unique()[1] - ground_truth_df['az_tru'].unique()[0]

fig, axs = plt.subplots(2, 1, figsize=(9, 10))

axs[0].bar(inc_bins[:-1], correct_proportions_inc, width=bar_width_inc, align='edge', label='Correct')
axs[0].bar(inc_bins[:-1], false_positives_proportions_inc, width=bar_width_inc, align='edge', bottom=correct_proportions_inc, label='False Positives')
axs[0].bar(inc_bins[:-1], false_negatives_proportions_inc, width=bar_width_inc, align='edge', bottom=correct_proportions_inc + false_positives_proportions_inc, label='False Negatives')
axs[0].set_xlabel('θ, °')
axs[0].set_ylabel('Proportion, %')
axs[0].set_title('θ')
axs[0].legend(loc='lower right')

axs[1].bar(az_bins[:-1], correct_proportions_az, width=bar_width_az, align='edge', label='Correct')
axs[1].bar(az_bins[:-1], false_positives_proportions_az, width=bar_width_az, align='edge', bottom=correct_proportions_az, label='False Positives')
axs[1].bar(az_bins[:-1], false_negatives_proportions_az, width=bar_width_az, align='edge', bottom=correct_proportions_az + false_positives_proportions_az, label='False Negatives')
axs[1].set_xlabel('φ, °')
axs[1].set_ylabel('Proportion')
axs[1].set_title('φ')
axs[1].legend(loc='lower right')

#plt.show()
plt.savefig(f"accuracy_step1deg", dpi=300)
plt.close()
