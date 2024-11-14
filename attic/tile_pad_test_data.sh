#!/bin/bash

# Input and output file names
input_file="data.txt"
output_file="data_tiled_4x4.txt"

#cp $input_file $output_file

{
    # Add 6 rows of zeros at the top
    for i in {1..6}; do printf '0 %.0s' {1..30}; echo; done

    # Read the original file, add padding, and format
    while read -r line; do
        printf '0 %.0s' {1..6};  # Add 6 columns of zeros on the left
        printf '%s' "$line";      # Print the original row
        printf ' 0%.0s' {1..12};  # Add 12 columns of zeros on the right
        echo;                    # Newline
    done < $input_file

    # Add 12 rows of zeros at the bottom
    for i in {1..12}; do printf '0 %.0s' {1..30}; echo; done
} > blah1

# tile to 4x4
paste blah1 blah1 blah1 blah1 > blah2
cat blah2 blah2 blah2 blah2 > $output_file

rm blah*
