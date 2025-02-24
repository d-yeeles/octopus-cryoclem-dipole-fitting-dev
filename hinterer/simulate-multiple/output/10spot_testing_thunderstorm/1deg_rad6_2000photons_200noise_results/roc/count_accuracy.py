import pandas as pd
import numpy as np

thunderstorm_df = pd.read_csv('../thunderstorm_results_gaussian.csv')
ground_truth_df = pd.read_csv('../ground_truth_gaussian.csv')

radius = 2000

def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

correct = []
false_positives = []
false_negatives = []

for frame in ground_truth_df['frame'].unique():
    gt_points = ground_truth_df[ground_truth_df['frame'] == frame].copy()
    ts_points = thunderstorm_df[thunderstorm_df['frame'] == frame].copy()

    matched_estimates = set()
    matched_gt = set()

    # Find nearest estimated point for each ground truth point
    for _, gt_row in gt_points.iterrows():
        distances = ts_points.apply(lambda row: calculate_distance(gt_row['x_tru'], gt_row['y_tru'], row['x [nm]'], row['y [nm]']), axis=1)
        nearby_estimates = distances[distances <= radius]

        # If there are points within the radius, first find the nearest one
        if len(nearby_estimates) > 0:
            closest_estimate_idx = nearby_estimates.idxmin()
            closest_estimate = ts_points.loc[closest_estimate_idx]

            # Mark that as 'correct'
            correct_row = closest_estimate.copy()
            correct_row['az_tru'] = gt_row['az_tru']
            correct_row['inc_tru'] = gt_row['inc_tru']
            correct.append(correct_row)
            matched_estimates.add(closest_estimate_idx)
            matched_gt.add(gt_row.name) # Save the corresponding frame (so we know inc/az that worked/fialed)

            # Any additional estimates within the radius are false positives
            for idx in nearby_estimates.index:
                if idx != closest_estimate_idx:
                    false_positives_row = ts_points.loc[idx].copy()
                    false_positives_row['az_tru'] = gt_row['az_tru']
                    false_positives_row['inc_tru'] = gt_row['inc_tru']
                    false_positives.append(false_positives_row)

        # If there are no points within radius, mark as false negative
        if len(nearby_estimates) == 0:
            false_negatives_row = gt_row.copy()
            false_negatives_row['az_tru'] = gt_row['az_tru']
            false_negatives_row['inc_tru'] = gt_row['inc_tru']
            false_negatives.append(false_negatives_row)

    # Any remaining estimates (i.e. outside radius) are false positives
    for idx, ts_row in ts_points.iterrows():
        if idx not in matched_estimates:
            false_positives_row = ts_row.copy()
            false_positives_row['az_tru'] = np.nan
            false_positives_row['inc_tru'] = np.nan
            false_positives.append(false_positives_row)

pd.DataFrame(correct).to_csv('results/TP_gaussian_radius2000.csv', index=False)
pd.DataFrame(false_positives).to_csv('results/FP_gaussian_radius2000.csv', index=False)
pd.DataFrame(false_negatives).to_csv('results/FN_gaussian_radius2000.csv', index=False)





thunderstorm_df = pd.read_csv('../thunderstorm_results_hinterer.csv')
ground_truth_df = pd.read_csv('../ground_truth_hinterer.csv')

correct = []
false_positives = []
false_negatives = []

for frame in ground_truth_df['frame'].unique():
    gt_points = ground_truth_df[ground_truth_df['frame'] == frame].copy()
    ts_points = thunderstorm_df[thunderstorm_df['frame'] == frame].copy()

    matched_estimates = set()
    matched_gt = set()

    # Find nearest estimated point for each ground truth point
    for _, gt_row in gt_points.iterrows():
        distances = ts_points.apply(lambda row: calculate_distance(gt_row['x_tru'], gt_row['y_tru'], row['x [nm]'], row['y [nm]']), axis=1)
        nearby_estimates = distances[distances <= radius]

        # If there are points within the radius, first find the nearest one
        if len(nearby_estimates) > 0:
            closest_estimate_idx = nearby_estimates.idxmin()
            closest_estimate = ts_points.loc[closest_estimate_idx]

            # Mark that as 'correct'
            correct_row = closest_estimate.copy()
            correct_row['az_tru'] = gt_row['az_tru']
            correct_row['inc_tru'] = gt_row['inc_tru']
            correct.append(correct_row)
            matched_estimates.add(closest_estimate_idx)
            matched_gt.add(gt_row.name) # Save the corresponding frame (so we know inc/az that worked/fialed)

            # Any additional estimates within the radius are false positives
            for idx in nearby_estimates.index:
                if idx != closest_estimate_idx:
                    false_positives_row = ts_points.loc[idx].copy()
                    false_positives_row['az_tru'] = gt_row['az_tru']
                    false_positives_row['inc_tru'] = gt_row['inc_tru']
                    false_positives.append(false_positives_row)

        # If there are no points within radius, mark as false negative
        if len(nearby_estimates) == 0:
            false_negatives_row = gt_row.copy()
            false_negatives_row['az_tru'] = gt_row['az_tru']
            false_negatives_row['inc_tru'] = gt_row['inc_tru']
            false_negatives.append(false_negatives_row)

    # Any remaining estimates (i.e. outside radius) are false positives
    for idx, ts_row in ts_points.iterrows():
        if idx not in matched_estimates:
            false_positives_row = ts_row.copy()
            false_positives_row['az_tru'] = np.nan
            false_positives_row['inc_tru'] = np.nan
            false_positives.append(false_positives_row)

pd.DataFrame(correct).to_csv('results/TP_hinterer_radius2000.csv', index=False)
pd.DataFrame(false_positives).to_csv('results/FP_hinterer_radius2000.csv', index=False)
pd.DataFrame(false_negatives).to_csv('results/FN_hinterer_radius2000.csv', index=False)

