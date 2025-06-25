#!/bin/bash

# Usage: ./update_toctree_refs.sh <directory>
directory="$1"

if [[ -z "$directory" ]]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Build a mapping of oldname-without-ext to newname-without-ext
declare -A name_map

for file in "$directory"/[0-9][0-9]-*.rst; do
    [[ -e "$file" ]] || continue
    basename=$(basename "$file" .rst)
    prefix="${basename:0:2}"
    rest="${basename:3}"
    name_map["$rest"]="$prefix-$rest"
done

# Update all .rst files
find "$directory" -type f -name '*.rst' | while read -r file; do
    tmpfile="${file}.tmp"

    cp "$file" "$tmpfile"

    for oldrest in "${!name_map[@]}"; do
        newbase="${name_map[$oldrest]}"
        # Replace only whole lines or indented lines containing the old base
        sed -i -E "s/^(\s*)$oldrest(\s*)$/\1$newbase\2/" "$tmpfile"
    done

    mv "$tmpfile" "$file"
done
