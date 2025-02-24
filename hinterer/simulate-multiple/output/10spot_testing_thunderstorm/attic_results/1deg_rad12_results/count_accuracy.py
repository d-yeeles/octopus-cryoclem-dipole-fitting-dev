import pandas as pd
import numpy as np

# Load the CSV files
thunderstorm_df = pd.read_csv('thunderstorm_results_gaussian.csv')
ground_truth_df = pd.read_csv('ground_truth_gaussian.csv')

# Define the radius for matching points (500 nm)
radius = 200

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
    # If there is no GT point within radius, add that GT point to 'false positives'

    for _, ts_row in ts_points.iterrows():
        distances = gt_points.apply(lambda row: calculate_distance(ts_row['x [nm]'], ts_row['y [nm]'], row['x_tru'], row['y_tru']), axis=1)

        if (distances <= radius).sum() == 1:
            correct.append(gt_points.loc[distances.idxmin()])  # Save the GT row corresponding to the nearest TS point
        elif (distances <= radius).sum() == 0:
            false_positives.append(gt_points.loc[distances.idxmin()])  # Save the GT row corresponding to the nearest TS point

    # Next, cycle through each GT point and compare it to every TS point
    # If there is more than 1 TS point within radius, add that GT point to 'false positives'
    # If there is no TS point within radius, add that GT point to 'false negatives'

    for _, gt_row in gt_points.iterrows():
        distances = ts_points.apply(lambda row: calculate_distance(gt_row['x_tru'], gt_row['y_tru'], row['x [nm]'], row['y [nm]']), axis=1)

        if (distances <= radius).sum() > 1:
            false_positives.append(gt_row)
        elif (distances <= radius).sum() == 0:
            false_negatives.append(gt_row)  # Save the GT row corresponding to the nearest TS point

pd.DataFrame(correct).to_csv('correct_gaussian.csv', index=False)
pd.DataFrame(false_positives).to_csv('false_positives_gaussian.csv', index=False)
pd.DataFrame(false_negatives).to_csv('false_negatives_gaussian.csv', index=False)

print("Processing complete.")












# Load the CSV files
thunderstorm_df = pd.read_csv('thunderstorm_results_hinterer.csv')
ground_truth_df = pd.read_csv('ground_truth_hinterer.csv')

# Define the radius for matching points (500 nm)
radius = 200

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
    # If there is no GT point within radius, add that GT point to 'false positives'

    for _, ts_row in ts_points.iterrows():
        distances = gt_points.apply(lambda row: calculate_distance(ts_row['x [nm]'], ts_row['y [nm]'], row['x_tru'], row['y_tru']), axis=1)

        if (distances <= radius).sum() == 1:
            correct.append(gt_points.loc[distances.idxmin()])  # Save the GT row corresponding to the nearest TS point
        elif (distances <= radius).sum() == 0:
            false_positives.append(gt_points.loc[distances.idxmin()])  # Save the GT row corresponding to the nearest TS point

    # Next, cycle through each GT point and compare it to every TS point
    # If there is more than 1 TS point within radius, add that GT point to 'false positives'
    # If there is no TS point within radius, add that GT point to 'false negatives'

    for _, gt_row in gt_points.iterrows():
        distances = ts_points.apply(lambda row: calculate_distance(gt_row['x_tru'], gt_row['y_tru'], row['x [nm]'], row['y [nm]']), axis=1)

        if (distances <= radius).sum() > 1:
            false_positives.append(gt_row)
        elif (distances <= radius).sum() == 0:
            #false_negatives.append(gt_points.loc[distances.idxmin()])  # Save the GT row corresponding to the nearest TS point
            false_negatives.append(gt_row)  # Save the GT row corresponding to the nearest TS point

pd.DataFrame(correct).to_csv('correct_hinterer.csv', index=False)
pd.DataFrame(false_positives).to_csv('false_positives_hinterer.csv', index=False)
pd.DataFrame(false_negatives).to_csv('false_negatives_hinterer.csv', index=False)

print("Processing complete.")





