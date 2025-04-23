#!/bin/bash

# Pad frame numbers in params_frame*.m files
for file in params_frame*.m; do
    [[ -e "$file" ]] || continue  # skip if no matching files
    num=$(echo "$file" | grep -oP '\d+')
    new=$(printf "params_frame%06d.m" "$num")
    mv "$file" "$new"
done

# Pad frame numbers in sim_frame*.tif files
for file in sim_frame*.tif; do
    [[ -e "$file" ]] || continue
    num=$(echo "$file" | grep -oP '\d+')
    new=$(printf "sim_frame%06d.tif" "$num")
    mv "$file" "$new"
done

