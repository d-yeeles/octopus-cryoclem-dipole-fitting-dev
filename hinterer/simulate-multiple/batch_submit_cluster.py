import csv
import subprocess
import argparse
import os
import OctopusCluster.SimpleCluster as Cluster

def get_frame_arrays(csv_path, chunk_size):
    # Increase the CSV field size limit
    csv.field_size_limit(1000000)
    unique_values = []
    seen = set()

    with open(csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header
        for row in reader:
            if len(row) >= 1 and row[0] not in seen:
                unique_values.append(row[0])
                seen.add(row[0])

    # Split into chunks of 100
    chunks = [unique_values[i:i+chunk_size] for i in range(0, len(unique_values), chunk_size)]

    return chunks

def main():
    parser = argparse.ArgumentParser(description='Run test thunderstorm script on frame chunks.')
    parser.add_argument('--image_path', required=True, help='Path to input image')
    parser.add_argument('--thunderstorm_path', required=True, help='Path to thunderstorm CSV')
    parser.add_argument('--output_path', required=True, help='Path to output results')
    parser.add_argument('--model', required=True, help='Model name')
    parser.add_argument('--patch_width', required=True, type=int, help='Patch width in nm')
    parser.add_argument('--chunk_size', required=True, type=int, help='How many frames to chunk the job up into')

    args = parser.parse_args()

    frame_array_chunks = get_frame_arrays(args.thunderstorm_path, args.chunk_size)

    # Extract base filename and directory from output path
    output_dir = os.path.dirname(args.output_path)
    output_filename = os.path.basename(args.output_path)
    filename_base, extension = os.path.splitext(output_filename)

    for i, frame_array in enumerate(frame_array_chunks):

        numeric_frames = [int(frame) for frame in frame_array]
        first_frame = min(numeric_frames)
        last_frame = max(numeric_frames)

        # Create a new output path with _chunkN appended
        chunk_output_path = os.path.join(output_dir, f"{filename_base}_chunk{i+1}{extension}")
        
        frame_string = ','.join(frame_array)  # Elements are already strings
        frame_arg = f'[{frame_string}]'       # Enclose in square brackets
        
        jobids=Cluster.submitCommands(
            Cluster.Cluster.octopus_cloud2,
            [[
                '/mnt/rclsfserv005/users/tfq96423/dipole_fitting/run_test_thunderstorm.sh',
                '/opt/matlabruntime/R2024b',
                args.image_path,
                args.thunderstorm_path,
                chunk_output_path,  # Use the new output path with chunk number
                args.model,
                str(args.patch_width),
                str(first_frame),
                str(last_frame),
            ]],
        
            singularity_image='/mnt/rclsfserv005/local/dipole_fitting.sif',
                 delete_successful_logs=True,
                 debug=True,
            )
        
        print(f'Submitted chunk {i+1} to {chunk_output_path}')

    print('All chunks submitted')

if __name__ == '__main__':
    main()
