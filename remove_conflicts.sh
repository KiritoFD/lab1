#!/bin/bash

echo "=== Removing Conflicting Files ==="

# Identify the C++ implementation of dna_repeat_finder
if [ -f "src/dna_repeat_finder.cpp" ]; then
    echo "Moving C++ implementation to cpp/ directory..."
    mkdir -p cpp
    mv src/dna_repeat_finder.cpp cpp/
    echo "  - Moved src/dna_repeat_finder.cpp to cpp/"
fi

# Handle conflicting dna_repeat_finder.c
if [ -f "src/dna_repeat_finder.c" ]; then
    echo "Checking for conflicts in dna_repeat_finder.c..."
    if grep -q "int main" src/dna_repeat_finder.c; then
        echo "  - Found main function in dna_repeat_finder.c"
        mkdir -p archive
        mv src/dna_repeat_finder.c archive/
        echo "  - Moved src/dna_repeat_finder.c to archive/"
    fi
fi

# Rename any duplicate header files in src/ when they exist in include/
for header in dna_common.h dna_io.h dna_traditional.h dna_graph.h; do
    if [ -f "src/$header" ] && [ -f "include/$header" ]; then
        echo "Found duplicate header: $header"
        mkdir -p archive
        mv src/$header archive/
        echo "  - Moved src/$header to archive/"
    fi
done

# Check for inconsistent or invalid include paths in source files
echo "Checking for include path issues..."
for src_file in src/*.c; do
    if grep -q "include \".*\.h\"" "$src_file"; then
        echo "  - Updating includes in $src_file"
        sed -i 's|include "dna_|include "../include/dna_|g' "$src_file"
    fi
done

echo ""
echo "Conflict resolution complete. Run the fix_project.sh script next to update the build system."
echo "=== Done ==="
