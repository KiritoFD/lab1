#!/bin/bash

echo "=== Fixing Directory Structure ==="

# Check if source files are actually in src/core/
if [ -d "src/core" ]; then
    echo "Found src/core directory, fixing include paths..."
    for src_file in src/core/*.c; do
        if [ -f "$src_file" ]; then
            echo "  - Updating include paths in $src_file"
            sed -i 's|../include/|../../include/|g' "$src_file"
        fi
    done
else
    # Create proper directory structure
    echo "Creating proper directory structure..."
    mkdir -p src/core src/dag
    
    # Move C files to src/core if they're in the main src directory
    for file in src/*.c; do
        if [[ "$file" != "src/main.c" && "$file" != "src/generate_test_sequences.c" ]]; then
            filename=$(basename "$file")
            if [ -f "$file" ]; then
                echo "  - Moving $file to src/core/$filename"
                mv "$file" "src/core/$filename"
            fi
        fi
    done
    
    # Update Makefile to use the right paths
    if [ -f "Makefile" ]; then
        echo "Updating Makefile paths..."
        sed -i 's|$(SRC_DIR)/dna_|$(SRC_DIR)/core/dna_|g' Makefile
        sed -i 's|$(SRC_DIR)/core/main.c|$(SRC_DIR)/main.c|g' Makefile
    fi
fi

# Fix main.c includes if needed
if [ -f "src/main.c" ]; then
    echo "Fixing includes in main.c..."
    sed -i 's|../include/core/|../include/|g' src/main.c
fi

echo "=== Done ==="
