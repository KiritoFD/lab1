#!/bin/bash

echo "=== DNA Repeat Finder Project Cleanup ==="

# Clean object files
echo "Cleaning object files..."
rm -rf obj/*.o obj/*.d

# Clean binaries
echo "Cleaning binaries..."
rm -rf bin/*

# Clean conflicts - handle the generate_test_sequences.c file
echo "Resolving file conflicts..."
if [ -f src/generate_test_sequences.c ]; then
    echo "  - Moving generate_test_sequences.c to tools/ directory"
    mkdir -p tools
    mv src/generate_test_sequences.c tools/
fi

# Move C++ implementation to a separate directory for clarity
if [ -f src/dna_repeat_finder.cpp ]; then
    echo "  - Moving C++ implementation to cpp/ directory"
    mkdir -p cpp
    mv src/dna_repeat_finder.cpp cpp/
fi

# Organize DAG-related files
if [ -f dag_algorithms.c ] || [ -f dag_algorithms.h ]; then
    echo "  - Organizing DAG algorithm files"
    mkdir -p algorithms
    [ -f dag_algorithms.c ] && mv dag_algorithms.c algorithms/
    [ -f dag_algorithms.h ] && mv dag_algorithms.h algorithms/
    [ -f dag_path_planning.md ] && mv dag_path_planning.md algorithms/
fi

# Clean temporary files
echo "Cleaning temporary and backup files..."
find . -name "*~" -delete
find . -name "*.bak" -delete
find . -name ".*.swp" -delete

echo "=== Cleanup Complete ==="
echo "You may now rebuild the project with: make clean && make"
