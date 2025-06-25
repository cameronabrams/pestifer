#!/bin/bash

# Usage: ./increment_filenames.sh <start_number>
start_from="$1"

if [[ -z "$start_from" ]]; then
    echo "Usage: $0 <start_number>"
    exit 1
fi

# Pad the starting number to two digits
printf -v start_padded "%02d" "$start_from"

# Gather all matching files with 2-digit numeric prefix, and sort in reverse order
files=($(ls [0-9][0-9]-* 2>/dev/null | sort -r))

for file in "${files[@]}"; do
    prefix="${file:0:2}"
    rest="${file:3}"

    # Compare with start_padded
    if (( 10#$prefix >= 10#$start_padded )); then
        new_index=$((10#$prefix + 1))
        printf -v new_prefix "%02d" "$new_index"
        new_name="${new_prefix}-${rest}"
        echo "Renaming: $file -> $new_name"
        mv "$file" "$new_name"
    fi
done
