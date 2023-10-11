#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 inputfile outputfile"
    exit 1
fi

awk -F'\t' 'BEGIN { count = 0; } { empty = 1; for (i = 9; i <= 19; i++) { if ($i != "") { empty = 0; break; } } if (empty == 1) { count++; } } END { print count; }' "$1" > "$2"
