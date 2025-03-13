#!/bin/bash

echo "=== Fixing Include Paths ==="

# Fix includes in all core source files
if [ -d "src/core" ]; then
  echo "Processing src/core files..."
  find src/core -name "*.c" -type f -exec sed -i 's|#include "../core/dna_common.h"|#include "../../include/dna_common.h"|g' {} \;
  find src/core -name "*.c" -type f -exec sed -i 's|#include "../include/dna_common.h"|#include "../../include/dna_common.h"|g' {} \;
  find src/core -name "*.c" -type f -exec sed -i 's|#include "../include/dna_io.h"|#include "../../include/dna_io.h"|g' {} \;
  find src/core -name "*.c" -type f -exec sed -i 's|#include "../include/dna_traditional.h"|#include "../../include/dna_traditional.h"|g' {} \;
  find src/core -name "*.c" -type f -exec sed -i 's|#include "../include/dna_graph.h"|#include "../../include/dna_graph.h"|g' {} \;
  find src/core -name "*.c" -type f -exec sed -i 's|#include "../include/cpu_optimize.h"|#include "../../include/cpu_optimize.h"|g' {} \;
fi

echo "=== Include Paths Fixed ==="
