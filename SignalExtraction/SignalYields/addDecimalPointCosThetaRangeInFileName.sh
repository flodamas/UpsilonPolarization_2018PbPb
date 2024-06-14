#!/bin/bash

# Directory containing the files
directory="/Users/ahram/Codes/Polarization/UpsilonPolarization_2018PbPb/SignalExtraction/SignalYields"

# Loop through each file in the directory
for file in "$directory"/*cosTheta*; do
    # Extract parts of the filename using pattern matching
    base="${file%_cosTheta*}"
    suffix="${file##*_cosTheta}"

    # Extract the cosTheta range
    cosTheta_range="${suffix%%_*}"

    # Check if the cosTheta range already has exactly two decimal points
    if [[ "$cosTheta_range" =~ \.[0-9]*\.[0-9]* ]]; then
        echo "Skipping $file: Already has two decimal points"
        continue
    fi


    # Extract start and end values
    start="${cosTheta_range%to*}"
    end="${cosTheta_range#*to}"

    # Format start and end values to two decimal points
    formatted_start=$(printf "%.2f" "${start//p/.}")
    formatted_end=$(printf "%.2f" "${end//p/.}")

    # Construct the new filename
    new_filename="${base}_cosTheta-${formatted_start}to${formatted_end}_${suffix#*_cosTheta*}"

    # Rename the file
    mv "$file" "$new_filename"
done
