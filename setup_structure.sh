#!/bin/bash

# Script to set up project directory structure
echo "Setting up DNA Repeat Finder project structure..."

# Create directories
mkdir -p src/core
mkdir -p src/tools
mkdir -p src/dag
mkdir -p include/core
mkdir -p include/dag
mkdir -p cpp
mkdir -p bin
mkdir -p obj
mkdir -p data
mkdir -p docs

# Move core C files
mv src/dna_common.c src/core/ 2>/dev/null
mv src/dna_io.c src/core/ 2>/dev/null
mv src/dna_traditional.c src/core/ 2>/dev/null
mv src/dna_graph.c src/core/ 2>/dev/null
mv src/dna_sequence_utils.c src/core/ 2>/dev/null

# Move tool files
mv src/generate_test_sequences.c src/tools/ 2>/dev/null

# Move DAG files
mv dag_algorithms.c src/dag/ 2>/dev/null
mv dag_path_planning.md docs/ 2>/dev/null

# Move header files
mv include/dna_common.h include/core/ 2>/dev/null
mv include/dna_io.h include/core/ 2>/dev/null
mv include/dna_traditional.h include/core/ 2>/dev/null
mv include/dna_graph.h include/core/ 2>/dev/null
mv include/cpu_optimize.h include/core/ 2>/dev/null
mv dag_algorithms.h include/dag/ 2>/dev/null

# Move C++ files
mv src/dna_repeat_finder.cpp cpp/ 2>/dev/null

# Set up data directory
touch data/reference.txt
touch data/query.txt

echo "Directory structure setup complete!"
