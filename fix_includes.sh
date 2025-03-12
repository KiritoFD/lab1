#!/bin/bash

echo "=== Fixing Include Paths ==="

# Fix include paths in src/core files
if [ -d "src/core" ]; then
    for src_file in src/core/*.c; do
        if [ -f "$src_file" ]; then
            filename=$(basename "$src_file")
            echo "Fixing includes in $filename..."
            sed -i 's|../include/|../../include/|g' "$src_file"
        fi
    done
fi

# Fix include paths in regular src files
for src_file in src/*.c; do
    if [ -f "$src_file" ]; then
        filename=$(basename "$src_file")
        echo "Fixing includes in $filename..."
        sed -i 's|../include/core/|../include/|g' "$src_file"
    fi
done

echo "=== Include Paths Fixed ==="
