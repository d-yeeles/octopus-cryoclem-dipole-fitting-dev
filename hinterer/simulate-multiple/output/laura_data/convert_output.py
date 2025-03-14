'''
USAGE: python laura_to_dave.py "string_match*" output.py
'''

import pandas as pd
import sys
import glob

def combine_csv_files(input_pattern):
    # Find all CSV files matching the pattern
    files = glob.glob(input_pattern)

    if not files:
        print(f"No files found matching the pattern: {input_pattern}")
        sys.exit(1)

    # Read the first file to use its header
    df_combined = pd.read_csv(files[0])

    # Loop through the rest of the files and append them
    for file in files[1:]:
        df = pd.read_csv(file)
        df_combined = pd.concat([df_combined, df], ignore_index=True)

    # Return the combined DataFrame
    return df_combined

def csv_to_python(df, output_py_file):
    # Initialize an empty string to store the formatted output
    python_output = ""

    # Iterate through each column and format it like the example
    for col in df.columns:
        # Replace words as specified
        new_col = col.replace('true', 'tru').replace('error', 'err').replace('theta', 'inc').replace('phi', 'az').replace('photons', 'photon')
        values = df[col].tolist()
        python_output += f"{new_col} = {values}\n"

    # Append additional lines
    python_output += "\nobj_tru = []\nobj_est = []\nobj_err = []\n"

    # Write to the output Python file
    with open(output_py_file, "w") as f:
        f.write(python_output)

    print(f"Converted the combined data to {output_py_file}")

# Check for correct number of arguments
if len(sys.argv) != 3:
    print("Usage: python blah.py input_pattern output.py")
    sys.exit(1)

# Get file names from command line arguments
input_pattern = sys.argv[1]
output_py = sys.argv[2]

# Combine the CSV files into one DataFrame
df_combined = combine_csv_files(input_pattern)

# Run the conversion on the combined DataFrame
csv_to_python(df_combined, output_py)

