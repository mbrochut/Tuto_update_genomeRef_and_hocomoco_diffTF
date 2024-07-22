#!/bin/bash

# Define the directory containing the .bed files
directory="../results/PWMScan_hocomocov12"

# Loop through each .bed file in the directory
for file in "$directory"/*.bed; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Use grep to filter lines that begin with chr[1-19] or chr[X,Y]
        grep -E '^chr([1-9]|1[0-9]|X|Y)\b' "$file" > "${file}.tmp"
        
        # Overwrite the original file with the filtered content
        mv "${file}.tmp" "$file"
        
        echo "Processed $file"
    fi
done