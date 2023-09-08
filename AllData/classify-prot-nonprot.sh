#!/bin/bash

# Create the "Proteins" folder if it doesn't exist
mkdir -p Proteins

# Iterate over all files in the current directory
for file in *; do
    # Check if the file is a regular file
    if [[ -f $file ]]; then
        # Check if the line contains "COMPND   2 MOLECULE:" and "PROTEIN"
        if grep -q "COMPND   2 MOLECULE:" "$file" && grep -q "PROTEIN" "$file"; then
            # Copy the file to the "Proteins" folder
            cp "$file" Proteins/
        fi
    fi
done

