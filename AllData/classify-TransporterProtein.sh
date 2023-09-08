#!/bin/bash

# Create the "Proteins" folder if it doesn't exist
foldername="Transporter-Proteins"
mkdir -p $foldername

# Iterate over all files in the current directory
for file in *; do
    # Check if the file is a regular file
    if [[ -f $file ]]; then
        # Check if the line contains "COMPND   2 MOLECULE:" and "PROTEIN"
        if grep -q "TRANSPORT PROTEIN" "$file"; then
            # Copy the file to the "Proteins" folder
            cp "$file" $foldername/
        fi
    fi
done

