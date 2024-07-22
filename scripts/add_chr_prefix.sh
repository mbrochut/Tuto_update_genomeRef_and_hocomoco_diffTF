#!/bin/bash

# Define the directory containing the files
directory="../results/PWMScan_hocomocov12"

# Loop through each file in the directory
for file in "$directory"/*; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Use sed to add "chr" in front of every line and overwrite the file
        sed -i 's/^/chr/' "$file"
        echo "Updated $file"
    fi
done

