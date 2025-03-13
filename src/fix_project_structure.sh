#!/bin/bash

echo "=== Fixing DNA Repeat Finder Project Structure ==="

# 1. Create necessary directories
mkdir -p src include obj bin

# 2. Delete main.c in root directory if it exists
if [ -f main.c ]; then
  echo "Removing main.c from root directory"
  rm main.c
fi

# 3. Move any stranded code files to src directory
if [ -f dna_common.c ]; then mv dna_common.c src/; fi
if [ -f dna_io.c ]; then mv dna_io.c src/; fi
if [ -f dna_traditional.c ]; then mv dna_traditional.c src/; fi
if [ -f dna_graph.c ]; then mv dna_graph.c src/; fi
if [ -f dna_sequence_utils.c ]; then mv dna_sequence_utils.c src/; fi
if [ -f dag_algorithms.c ]; then mv dag_algorithms.c src/; fi

# 4. Move any stranded header files to include directory
if [ -f dna_common.h ]; then mv dna_common.h include/; fi
if [ -f dna_io.h ]; then mv dna_io.h include/; fi
if [ -f dna_traditional.h ]; then mv dna_traditional.h include/; fi
if [ -f dna_graph.h ]; then mv dna_graph.h include/; fi
if [ -f cpu_optimize.h ]; then mv cpu_optimize.h include/; fi
if [ -f dag_algorithms.h ]; then mv dag_algorithms.h include/; fi

# 5. Fix any troublesome include paths in src directory files
find src -name "*.c" -type f -exec sed -i 's|"../include/|"../include/|g' {} \;
find src -name "*.c" -type f -exec sed -i 's|"../core/dna_|"../include/dna_|g' {} \;
find src -name "*.c" -type f -exec sed -i 's|"dna_common.h"|"../include/dna_common.h"|g' {} \;
find src -name "*.c" -type f -exec sed -i 's|"dna_io.h"|"../include/dna_io.h"|g' {} \;
find src -name "*.c" -type f -exec sed -i 's|"dna_traditional.h"|"../include/dna_traditional.h"|g' {} \;
find src -name "*.c" -type f -exec sed -i 's|"dna_graph.h"|"../include/dna_graph.h"|g' {} \;
find src -name "*.c" -type f -exec sed -i 's|"cpu_optimize.h"|"../include/cpu_optimize.h"|g' {} \;

# 6. Fix any include paths in src/core directory if it exists
if [ -d "src/core" ]; then
  find src/core -name "*.c" -type f -exec sed -i 's|"../include/|"../../include/|g' {} \;
  find src/core -name "*.c" -type f -exec sed -i 's|"../core/dna_|"../../include/dna_|g' {} \;
fi

# Special case for src/core/dna_common.c
if [ -f "src/core/dna_common.c" ]; then
  echo "Fixing src/core/dna_common.c include paths"
  sed -i '1s|#include "../include/dna_common.h"|#include "../../include/dna_common.h"|' src/core/dna_common.c
  sed -i 's|#include "../core/dna_common.h"|#include "../../include/dna_common.h"|g' src/core/dna_common.c
fi

echo "=== Project Structure Fixed ==="
echo "Now run: make clean && make"
