#!/bin/bash

# Define the directory containing the .bed files
directory="../results/PWMScan_hocomocov12"

# Loop through all .bed files in the directory
for file in "$directory"/*.bed; do
    # Extract the base name without .bed
    base_name=$(basename ${file})
    # Use cut to get the NAME part
    name=$(echo "$base_name" | cut -d '.' -f 1)
    # Construct the new file name
    new_name="${directory}/${name}_TFBS.bed"
    # Rename the file
    mv "$file" "$new_name"
done