import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os
from collections import defaultdict

def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

def compute_dice_for_frames_by_inclination(thunderstorm_csv, groundtruth_csv, inclination_range, radius=200):
    """
    Compute DICE coefficient for frames with molecules in the specified inclination range.
    
    Parameters:
    thunderstorm_csv: path to the ThunderSTORM results CSV
    groundtruth_csv: path to the ground truth CSV
    inclination_range: tuple of (min_inc, max_inc) for the inclination range in degrees
    radius: maximum matching distance in nm
    
    Returns:
    Dictionary of DICE results for each frame in the specified inclination range
    """
    # Load data
    thunderstorm_df = pd.read_csv(thunderstorm_csv)
    ground_truth_df = pd.read_csv(groundtruth_csv)
    
    # Convert inclination from radians to degrees
    if 'inc_tru' in ground_truth_df.columns:
        ground_truth_df['inc_tru_deg'] = ground_truth_df['inc_tru'] * (180 / np.pi)
    else:
        print("Warning: 'inc_tru' column not found in ground truth data")
        return {}
    
    # Filter ground truth data by inclination range
    min_inc, max_inc = inclination_range
    ground_truth_filtered = ground_truth_df[
        (ground_truth_df['inc_tru_deg'] >= min_inc) & 
        (ground_truth_df['inc_tru_deg'] < max_inc)
    ]
    
    # Get unique frames containing molecules in this inclination range
    frames_in_range = ground_truth_filtered['frame'].unique()
    
    # Initialize results
    dice_results = {}
    
    # Process each frame
    for frame in frames_in_range:
        # Filter ground truth points for this frame and inclination range
        gt_points = ground_truth_filtered[ground_truth_filtered['frame'] == frame].copy()
        
        # Get ThunderSTORM points for this frame
        ts_points = thunderstorm_df[thunderstorm_df['frame'] == frame].copy()
        
        # Skip if either dataset is empty for this frame
        if len(gt_points) == 0 or len(ts_points) == 0:
            continue
            
        matched_estimates = set()
        
        tp = 0  # True positives
        fn = 0  # False negatives
        
        # Find nearest estimated point for each ground truth point
        for _, gt_row in gt_points.iterrows():
            distances = ts_points.apply(
                lambda row: calculate_distance(
                    gt_row['x_tru'], gt_row['y_tru'], 
                    row['x [nm]'], row['y [nm]']
                ), axis=1
            )
            nearby_estimates = distances[distances <= radius]
            
            # If there are points within the radius, find the nearest one
            if len(nearby_estimates) > 0:
                closest_estimate_idx = nearby_estimates.idxmin()
                
                # Mark as true positive
                tp += 1
                matched_estimates.add(closest_estimate_idx)
            else:
                # No points within radius - false negative
                fn += 1
        
        # Any remaining estimates are false positives
        fp = len(ts_points) - len(matched_estimates)
        
        # Calculate Dice coefficient
        if tp == 0:
            dice = 0
        else:
            dice = (2 * tp) / (2 * tp + fp + fn)
        
        # Store results
        dice_results[frame] = {
            'dice': dice,
            'tp': tp,
            'fp': fp,
            'fn': fn
        }
    
    return dice_results

def compute_overall_dice(dice_results):
    """
    Compute the overall DICE coefficient for the entire dataset.
    """
    # Initialize counters
    total_tp = 0
    total_fp = 0
    total_fn = 0
    
    # Sum up all TP, FP, FN across frames
    for frame_data in dice_results.values():
        total_tp += frame_data['tp']
        total_fp += frame_data['fp']
        total_fn += frame_data['fn']
    
    # Calculate overall DICE coefficient
    if total_tp == 0:
        overall_dice = 0
    else:
        overall_dice = (2 * total_tp) / (2 * total_tp + total_fp + total_fn)
    
    # Return results
    overall_results = {
        'overall_dice': overall_dice,
        'total_tp': total_tp,
        'total_fp': total_fp,
        'total_fn': total_fn,
        'total_detections': total_tp + total_fp,
        'total_ground_truth': total_tp + total_fn
    }
    
    return overall_results

def parse_filename(filename):
    """
    Parse filename to extract model, sigma, and fitrad values.
    Format: "thunderstorm-[model]-sigma[valueA]-fitrad[valueB].csv"
    """
    # Extract base filename if it's a path
    base_filename = os.path.basename(filename)
    pattern = r'thunderstorm-([^-]+)-sigma([0-9.]+)-fitrad([0-9.]+)\.csv'
    match = re.match(pattern, base_filename)
    
    if match:
        model = match.group(1)
        sigma = float(match.group(2))
        fitrad = float(match.group(3))
        return model, sigma, fitrad
    else:
        return None, None, None

def find_groundtruth_file(model, groundtruth_files):
    """
    Find the corresponding ground truth file for a given model.
    """
    for gt_file in groundtruth_files:
        if model in gt_file:
            return gt_file
    
    # Default fallback if specific model ground truth not found
    return None

def process_all_files_by_inclination(thunderstorm_files, groundtruth_files, inclination_ranges, radius=200):
    """
    Process all thunderstorm result files and compute DICE coefficients for each inclination range.
    
    Parameters:
    thunderstorm_files: List of ThunderSTORM result CSV files
    groundtruth_files: List of ground truth CSV files
    inclination_ranges: List of (min_inc, max_inc) tuples defining the inclination ranges
    radius: Maximum matching distance in nm
    
    Returns:
    Dictionary mapping inclination ranges to lists of result dictionaries
    """
    # Initialize results dictionary
    results_by_inclination = {inc_range: [] for inc_range in inclination_ranges}
    
    for ts_file in thunderstorm_files:
        # Extract parameters from filename
        model, sigma, fitrad = parse_filename(ts_file)
        
        if model is None:
            print(f"Skipping file with unknown format: {ts_file}")
            continue
        
        # Find corresponding ground truth file
        gt_file = find_groundtruth_file(model, groundtruth_files)
        if gt_file is None:
            print(f"No ground truth file found for model {model}, skipping {ts_file}")
            continue
            
        print(f"Processing {ts_file} with {gt_file}")
        
        # Process each inclination range
        for inc_range in inclination_ranges:
            min_inc, max_inc = inc_range
            print(f"  Analyzing inclination range {min_inc}° to {max_inc}°")
            
            # Compute DICE results for this inclination range
            dice_results = compute_dice_for_frames_by_inclination(ts_file, gt_file, inc_range, radius)
            
            # Skip if no results for this range
            if not dice_results:
                print(f"  No frames found for inclination range {min_inc}° to {max_inc}° in {ts_file}")
                continue
                
            overall_results = compute_overall_dice(dice_results)
            
            # Store results with parameters
            result_entry = {
                'file': ts_file,
                'model': model,
                'sigma': sigma,
                'fitrad': fitrad,
                'inc_min': min_inc,
                'inc_max': max_inc,
                'dice': overall_results['overall_dice'],
                'tp': overall_results['total_tp'],
                'fp': overall_results['total_fp'],
                'fn': overall_results['total_fn']
            }
            
            results_by_inclination[inc_range].append(result_entry)
            
            print(f"    {model} (sigma={sigma}, fitrad={fitrad}): DICE = {overall_results['overall_dice']:.4f}")
            print(f"    TP={overall_results['total_tp']}, FP={overall_results['total_fp']}, FN={overall_results['total_fn']}")
    
    return results_by_inclination

def plot_heatmap_by_inclination(results, model, inc_range, output_filename):
    """
    Create a heatmap of DICE coefficients for different sigma and fitrad values
    for a specific model and inclination range.
    """
    # Get min and max inclination for display
    min_inc, max_inc = inc_range
    
    # Filter results for the given model
    model_results = [r for r in results if r['model'] == model]
    
    if not model_results:
        print(f"No results found for model {model} in inclination range {min_inc}° to {max_inc}°")
        return
    
    # Extract unique sigma and fitrad values
    sigma_values = sorted(set(r['sigma'] for r in model_results))
    fitrad_values = sorted(set(r['fitrad'] for r in model_results))
    
    # Check if we have complete data
    is_complete = True
    for sigma in sigma_values:
        for fitrad in fitrad_values:
            if not any(r['sigma'] == sigma and r['fitrad'] == fitrad for r in model_results):
                print(f"Warning: Missing data point for {model} (sigma={sigma}, fitrad={fitrad})")
                is_complete = False
    
    # Create empty matrix for heatmap
    heatmap_data = np.zeros((len(sigma_values), len(fitrad_values)))
    
    # Map values to indices
    sigma_to_idx = {val: idx for idx, val in enumerate(sigma_values)}
    fitrad_to_idx = {val: idx for idx, val in enumerate(fitrad_values)}
    
    # Fill matrix with DICE coefficients
    for r in model_results:
        i = sigma_to_idx[r['sigma']]
        j = fitrad_to_idx[r['fitrad']]
        heatmap_data[i, j] = r['dice']
    
    # Create figure with specified aspect ratio to make cells square
    plt.figure(figsize=(6, 8))
    
    # Plot the heatmap
    plt.imshow(heatmap_data, cmap='viridis', interpolation='nearest')
    
    # Add colorbar and labels
    cbar = plt.colorbar(label='Dice Coefficient')
    cbar.ax.tick_params(labelsize=9)
    
    plt.xlabel('Fitting Radius, px')
    plt.ylabel('Gaussian Blur Sigma, px')
    
    # Set ticks to parameter values
    plt.xticks(range(len(fitrad_values)), [str(v) for v in fitrad_values], rotation=45)
    plt.yticks(range(len(sigma_values)), [str(v) for v in sigma_values])
    
    # Find the best value
    max_idx = np.unravel_index(np.argmax(heatmap_data), heatmap_data.shape)
    
    # Add text annotations with the values
    for i in range(len(sigma_values)):
        for j in range(len(fitrad_values)):
            plt.text(j, i, f"{heatmap_data[i, j]:.3f}", 
                    ha="center", va="center", 
                    color="white" if heatmap_data[i, j] < 0.7 else "black",
                    fontsize=9)
    
    # Show the best parameter combination
    max_sigma = sigma_values[max_idx[0]]
    max_fitrad = fitrad_values[max_idx[1]]
    max_dice = heatmap_data[max_idx]
    
    plt.title(f'Dice Coefficient - {model.capitalize()}\nInclination {min_inc}° to {max_inc}°\n' r'$\frac{2TP}{2TP + FP + FN}$')
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    plt.close()

def check_for_new_files(existing_results):
    """
    Check if there are any new files that need to be processed
    compared to existing results
    """
    if not os.path.exists('./results'):
        print("Error: './results' directory not found!")
        return False, []
        
    result_files = os.listdir('./results')
    thunderstorm_files = [os.path.join('./results', f) for f in result_files if f.startswith('thunderstorm-') and 'sigma' in f and 'fitrad' in f]
    
    # Create a set of files already processed
    processed_files = set()
    for inc_results in existing_results.values():
        processed_files.update(set(r['file'] for r in inc_results))
    
    # Check for new files
    new_files = [f for f in thunderstorm_files if f not in processed_files]
    
    if new_files:
        print(f"Found {len(new_files)} new thunderstorm result files that need processing")
        return True, thunderstorm_files
    else:
        print("No new files found that need processing")
        return False, []

def main():
    # Define inclination ranges in degrees
    inclination_ranges = [
        (0, 22.5),
        (22.5, 45),
        (45, 67.5),
        (67.5, 90)
    ]
    
    print("Analyzing inclination ranges (in degrees):")
    for min_inc, max_inc in inclination_ranges:
        print(f"  {min_inc}° to {max_inc}°")
    
    # Variable to store results, either from file or newly computed
    results_by_inclination = {inc_range: [] for inc_range in inclination_ranges}
    need_computation = False
    thunderstorm_files = []
    
    # Check if we already have computed dice coefficients
    if os.path.exists('dc_rad200_by_inclination.csv'):
        print("Found existing dc_rad200_by_inclination.csv file. Loading saved results...")
        try:
            results_df = pd.read_csv('dc_rad200_by_inclination.csv')
            
            # Group results by inclination range
            for inc_range in inclination_ranges:
                min_inc, max_inc = inc_range
                range_results = results_df[
                    (results_df['inc_min'] == min_inc) & 
                    (results_df['inc_max'] == max_inc)
                ]
                results_by_inclination[inc_range] = range_results.to_dict('records')
                print(f"Loaded {len(results_by_inclination[inc_range])} pre-computed results for inclination range {min_inc}° to {max_inc}°")
            
            # Check if there are new files to process
            has_new_files, thunderstorm_files = check_for_new_files(results_by_inclination)
            
            if has_new_files:
                print("New files detected that aren't in dc_rad200_by_inclination.csv.")
                # Ask if user wants to compute results for new files
                recompute = input("Do you want to compute Dice coefficients for the new files? (y/n): ").lower()
                if recompute == 'y':
                    need_computation = True
                else:
                    print("Using only the saved results for plotting...")
            else:
                print("All files have been processed. Using saved results for plotting...")
        except Exception as e:
            print(f"Error loading dc_rad200_by_inclination.csv: {e}")
            print("Will recompute all results")
            need_computation = True
            
    else:
        print("No saved results found. Computing Dice coefficients...")
        need_computation = True
    
    # Compute new results if needed
    if need_computation:
        # If we don't have thunderstorm_files yet, get them
        if not thunderstorm_files:
            if not os.path.exists('./results'):
                print("Error: './results' directory not found!")
                return
                
            result_files = os.listdir('./results')
            thunderstorm_files = [os.path.join('./results', f) for f in result_files if f.startswith('thunderstorm-') and 'sigma' in f and 'fitrad' in f]
        
        # Find ground truth files in the current directory
        current_dir_files = os.listdir('.')
        groundtruth_files = [f for f in current_dir_files if f.startswith('ground_truth_')]
        
        print(f"Found {len(thunderstorm_files)} thunderstorm result files to process")
        print(f"Found {len(groundtruth_files)} ground truth files")
        
        # Process new files
        new_results_by_inclination = process_all_files_by_inclination(thunderstorm_files, groundtruth_files, inclination_ranges)
        
        # Merge new results with existing ones
        for inc_range, new_results in new_results_by_inclination.items():
            if not new_results:
                continue
                
            # If we already had results for this range, check for duplicates and update
            if results_by_inclination[inc_range]:
                # Create a dictionary to easily look up existing results by file, sigma, and fitrad
                existing_results = {
                    (r['file'], r['sigma'], r['fitrad']): i 
                    for i, r in enumerate(results_by_inclination[inc_range])
                }
                
                # Update existing results or add new ones
                for new_result in new_results:
                    key = (new_result['file'], new_result['sigma'], new_result['fitrad'])
                    if key in existing_results:
                        # Replace the existing result with the new one
                        results_by_inclination[inc_range][existing_results[key]] = new_result
                    else:
                        # Add the new result
                        results_by_inclination[inc_range].append(new_result)
            else:
                # If we didn't have any results for this range before, just use the new ones
                results_by_inclination[inc_range] = new_results
        
        # Save all results to CSV
        # First, flatten the results dictionary into a single list
        all_results = []
        for inc_range, range_results in results_by_inclination.items():
            all_results.extend(range_results)
        
        # Convert to DataFrame and save
        if all_results:
            results_df = pd.DataFrame(all_results)
            results_df.to_csv('dc_rad200_by_inclination.csv', index=False)
            print(f"Saved {len(all_results)} results to dc_rad200_by_inclination.csv")
    
    # Check if we have any results to plot
    has_results = False
    for inc_range, range_results in results_by_inclination.items():
        if range_results:
            has_results = True
            break
    
    if not has_results:
        print("No results to plot. Check file naming conventions and existence.")
        return
    
    # Plot with the results (either newly computed or loaded)
    plot_results(results_by_inclination, inclination_ranges)

def plot_results(results_by_inclination, inclination_ranges):
    """Generate all plots from the given results for each inclination range"""
    
    # Create a directory for the plots
    if not os.path.exists('./plots_by_inclination'):
        os.makedirs('./plots_by_inclination')
    
    # Create heatmaps for each model and inclination range
    for inc_range in inclination_ranges:
        min_inc, max_inc = inc_range
        range_results = results_by_inclination[inc_range]
        
        if not range_results:
            print(f"No results available for inclination range {min_inc}° to {max_inc}°. Skipping plots.")
            continue
        
        # Get unique models
        models = set(r['model'] for r in range_results)
        
        for model in models:
            # Filter results for this model and inclination range
            model_results = [r for r in range_results if r['model'] == model]
            
            # Generate heatmap
            output_filename = f'./plots_by_inclination/dice_heatmap_rad200_{model}_inc_{min_inc}to{max_inc}.png'
            plot_heatmap_by_inclination(model_results, model, inc_range, output_filename)
            print(f"Generated heatmap for {model}, inclination range {min_inc}° to {max_inc}°")
    
    print("All plots have been generated in the 'plots_by_inclination' directory.")

def plot_from_file():
    """
    Function to just generate plots from the saved dc_rad200_by_inclination.csv file
    without recomputing dice coefficients
    """
    if not os.path.exists('dc_rad200_by_inclination.csv'):
        print("Error: dc_rad200_by_inclination.csv not found! Please run the main function first to compute results.")
        return
    
    # Define inclination ranges in degrees
    inclination_ranges = [
        (0, 22.5),
        (22.5, 45),
        (45, 67.5),
        (67.5, 90)
    ]
    
    print("Analyzing inclination ranges (in degrees):")
    for min_inc, max_inc in inclination_ranges:
        print(f"  {min_inc}° to {max_inc}°")
    
    print("Loading results from dc_rad200_by_inclination.csv...")
    results_df = pd.read_csv('dc_rad200_by_inclination.csv')
    
    # Group results by inclination range
    results_by_inclination = {}
    for inc_range in inclination_ranges:
        min_inc, max_inc = inc_range
        range_results = results_df[
            (results_df['inc_min'] == min_inc) & 
            (results_df['inc_max'] == max_inc)
        ]
        results_by_inclination[inc_range] = range_results.to_dict('records')
        print(f"Loaded {len(results_by_inclination[inc_range])} results for inclination range {min_inc}° to {max_inc}°")
    
    plot_results(results_by_inclination, inclination_ranges)

if __name__ == "__main__":
    # Check if we're called with command line arguments
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--plot-only":
        print("Running in plot-only mode...")
        plot_from_file()
    else:
        main()
