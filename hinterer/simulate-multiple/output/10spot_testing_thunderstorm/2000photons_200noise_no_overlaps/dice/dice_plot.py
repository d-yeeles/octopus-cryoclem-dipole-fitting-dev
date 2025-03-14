import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os
from collections import defaultdict

def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

def compute_dice_for_frames(thunderstorm_csv, groundtruth_csv, radius=200):
    # Load data
    thunderstorm_df = pd.read_csv(thunderstorm_csv)
    ground_truth_df = pd.read_csv(groundtruth_csv)
    
    # Initialize results
    dice_results = {}
    
    # Process each frame
    for frame in ground_truth_df['frame'].unique():
        gt_points = ground_truth_df[ground_truth_df['frame'] == frame].copy()
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

def process_all_files(thunderstorm_files, groundtruth_files, radius=200):
    """
    Process all thunderstorm result files and compute DICE coefficients.
    """
    results = []
    
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
        
        # Compute DICE results
        dice_results = compute_dice_for_frames(ts_file, gt_file, radius)
        overall_results = compute_overall_dice(dice_results)
        
        # Store results with parameters
        result_entry = {
            'file': ts_file,
            'model': model,
            'sigma': sigma,
            'fitrad': fitrad,
            'dice': overall_results['overall_dice'],
            'tp': overall_results['total_tp'],
            'fp': overall_results['total_fp'],
            'fn': overall_results['total_fn']
        }
        
        results.append(result_entry)
        
        print(f"  {model} (sigma={sigma}, fitrad={fitrad}): DICE = {overall_results['overall_dice']:.4f}")
    
    return results

def plot_dice_vs_parameter(results, param_name, output_filename):
    """
    Plot DICE coefficient against a parameter (sigma or fitrad),
    grouping by model and other parameter values.
    """
    # Group results by model and the other parameter
    other_param = 'fitrad' if param_name == 'sigma' else 'sigma'
    grouped_results = defaultdict(list)
    
    for result in results:
        key = (result['model'], result[other_param])
        grouped_results[key].append((result[param_name], result['dice']))
    
    plt.figure(figsize=(12, 8))
    
    for (model, other_val), data_points in grouped_results.items():
        # Sort data points by parameter value
        data_points.sort(key=lambda x: x[0])
        
        # Extract sorted parameter values and DICE scores
        param_values, dice_scores = zip(*data_points)
        
        label = f"{model}, {other_param}={other_val}"
        plt.plot(param_values, dice_scores, marker='o', label=label, linewidth=2)
    
    plt.xlabel(f"{param_name} value")
    plt.ylabel("DICE Coefficient")
    plt.title(f"DICE Coefficient vs {param_name}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    plt.close()

def plot_heatmap(results, model, output_filename):
    """
    Create a heatmap of DICE coefficients for different sigma and fitrad values
    for a specific model.
    """
    # Filter results for the given model
    model_results = [r for r in results if r['model'] == model]
    
    if not model_results:
        print(f"No results found for model {model}")
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
    
    # Find the best and worst values
    max_idx = np.unravel_index(np.argmax(heatmap_data), heatmap_data.shape)
    min_idx = np.unravel_index(np.argmin(heatmap_data), heatmap_data.shape)
    
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
    
    plt.title(f'Dice Coefficient - {model.capitalize()}, 8100 spots\n' r'$\frac{2TP}{2TP + FP + FN}$')
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
    processed_files = set(r['file'] for r in existing_results) if existing_results else set()
    
    # Check for new files
    new_files = [f for f in thunderstorm_files if f not in processed_files]
    
    if new_files:
        print(f"Found {len(new_files)} new thunderstorm result files that need processing")
        return True, thunderstorm_files
    else:
        print("No new files found that need processing")
        return False, []

def main():
    # Variable to store results, either from file or newly computed
    results = []
    need_computation = False
    thunderstorm_files = []
    
    # Check if we already have computed dice coefficients
    if os.path.exists('dc_rad200.csv'):
        print("Found existing dc_rad200.csv file. Loading saved results...")
        try:
            results_df = pd.read_csv('dc_rad200.csv')
            results = results_df.to_dict('records')
            print(f"Loaded {len(results)} pre-computed results")
            
            # Check if there are new files to process
            has_new_files, thunderstorm_files = check_for_new_files(results)
            
            if has_new_files:
                print("New files detected that aren't in dc_rad200.csv.")
                # Ask if user wants to compute results for new files
                recompute = input("Do you want to compute Dice coefficients for the new files? (y/n): ").lower()
                if recompute == 'y':
                    need_computation = True
                else:
                    print("Using only the saved results for plotting...")
            else:
                print("All files have been processed. Using saved results for plotting...")
        except Exception as e:
            print(f"Error loading dc_rad200.csv: {e}")
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
        new_results = process_all_files(thunderstorm_files, groundtruth_files)
        
        if new_results:
            # If we already had results, check for duplicates and update
            if results:
                # Create a dictionary to easily look up existing results by file
                existing_files = {r['file']: i for i, r in enumerate(results)}
                
                # Update existing results or add new ones
                for new_result in new_results:
                    if new_result['file'] in existing_files:
                        # Replace the existing result with the new one
                        results[existing_files[new_result['file']]] = new_result
                    else:
                        # Add the new result
                        results.append(new_result)
            else:
                # If we didn't have any results before, just use the new ones
                results = new_results
            
            # Save all results to CSV
            results_df = pd.DataFrame(results)
            results_df.to_csv('dc_rad200.csv', index=False)
            print(f"Saved {len(results)} results to dc_rad200.csv")
    
    if not results:
        print("No results to plot. Check file naming conventions and existence.")
        return
    
    # Plot with the results (either newly computed or loaded)
    plot_results(results)

def plot_results(results):
    """Generate all plots from the given results"""
    
#    # Plot DICE vs sigma and DICE vs fitrad
#    plot_dice_vs_parameter(results, 'sigma', 'dice_vs_sigma.png')
#    plot_dice_vs_parameter(results, 'fitrad', 'dice_vs_fitrad.png')
    
    # Create heatmaps for each model
    models = set(r['model'] for r in results)
    for model in models:
        plot_heatmap(results, model, f'dice_heatmap_rad200_{model}.png')
    
    print("All plots have been generated.")

def plot_from_file():
    """
    Function to just generate plots from the saved dc_rad200.csv file
    without recomputing dice coefficients
    """
    if not os.path.exists('dc_rad200.csv'):
        print("Error: dc_rad200.csv not found! Please run the main function first to compute results.")
        return
    
    print("Loading results from dc_rad200.csv...")
    results_df = pd.read_csv('dc_rad200.csv')
    results = results_df.to_dict('records')
    print(f"Loaded {len(results)} results")
    
    plot_results(results)

if __name__ == "__main__":
    # Check if we're called with command line arguments
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--plot-only":
        print("Running in plot-only mode...")
        plot_from_file()
    else:
        main()
