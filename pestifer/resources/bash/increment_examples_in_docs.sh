#!/bin/bash
# written by chatgpt 4o
# Usage: ./increment_example_refs.sh <start_number> <directory>
start_from="$1"
directory="$2"

if [[ -z "$start_from" || -z "$directory" ]]; then
    echo "Usage: $0 <start_number> <directory>"
    exit 1
fi

# Walk the tree and process *.rst files
find "$directory" -type f -name '*.rst' | while read -r file; do
    echo "Processing $file"

    # Update both labels and refs
    perl -i -pe '
        use integer;
        s/(\.\. _example )(\d+)(:)/$1.($2 >= '"$start_from"' ? $2+1 : $2).$3/ge;
        s/(:ref:`example )(\d+)(`)/$1.($2 >= '"$start_from"' ? $2+1 : $2).$3/ge;
    ' "$file"
done

