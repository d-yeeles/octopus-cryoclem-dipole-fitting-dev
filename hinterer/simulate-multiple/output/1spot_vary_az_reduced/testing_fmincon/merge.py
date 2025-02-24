import re

def parse_file(filename):
    """ Reads a file and extracts variable names with their lists. """
    data = {}
    with open(filename, 'r') as file:
        for line in file:
            match = re.match(r'(\w+)\s*=\s*\[(.*)\]', line.strip())
            if match:
                var_name, values = match.groups()
                # Filter out empty strings and strip spaces before converting to float
                data[var_name] = [float(v.strip()) for v in values.split(',') if v.strip()]
    return data

def merge_data(file1, file2, output_file):
    """ Merges lists from two files and writes the combined output. """
    data1 = parse_file(file1)
    data2 = parse_file(file2)

    merged_data = {}
    for key in data1.keys():
        if key in data2:
            merged_data[key] = data1[key] + data2[key]
        else:
            merged_data[key] = data1[key]  # If a key exists in one file only

    with open(output_file, 'w') as f:
        for key, values in merged_data.items():
            values_str = ', '.join(f'{v:.4f}' for v in values)  # Format with 4 decimal places
            f.write(f"{key} = [{values_str} ]\n")  # Removed trailing comma

# Example usage
merge_data('fitting_results_hinterer_fmincon_180fudge_G.py-part1', 'fitting_results_hinterer_fmincon_180fudge_G.py-part2', 'fitting_results_hinterer_fmincon_180fudge_G.py')

