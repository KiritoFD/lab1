#!/bin/bash

echo "=== Creating Project Directory Structure ==="

# Clean up before creating the right structure
echo "Cleaning up before restructuring..."
rm -rf obj/*.o obj/*.d

# Create all needed directories
mkdir -p bin obj include src/core

# Check if we need to move files from src to src/core
if [ ! -d "src/core" ] || [ ! "$(ls -A src/core 2>/dev/null)" ]; then
    echo "Moving source files to src/core directory..."
    for file in src/dna_*.c; do
        if [ -f "$file" ]; then
            filename=$(basename "$file")
            echo "  - Moving $file to src/core/$filename"
            cp "$file" "src/core/$filename"
        fi
    done
fi

# Fix Makefile
echo "Updating Makefile..."
if [ -f "Makefile" ]; then
    # Create backup
    cp Makefile Makefile.backup
    
    # Replace source definitions to match structure
    sed -i 's|COMMON_SOURCES = $(SRC_DIR)/dna_common.c|COMMON_SOURCES = $(SRC_DIR)/core/dna_common.c|' Makefile
    sed -i 's|$(SRC_DIR)/dna_io.c|$(SRC_DIR)/core/dna_io.c|' Makefile
    sed -i 's|$(SRC_DIR)/dna_traditional.c|$(SRC_DIR)/core/dna_traditional.c|' Makefile
    sed -i 's|$(SRC_DIR)/dna_graph.c|$(SRC_DIR)/core/dna_graph.c|' Makefile
    sed -i 's|$(SRC_DIR)/dna_sequence_utils.c|$(SRC_DIR)/core/dna_sequence_utils.c|' Makefile
    
    # Fix compile rule for core files
    sed -i '/^# Compile core source files/,+1c\
# Compile core source files\
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c\
\t$(CC) $(CFLAGS) -MMD -c $< -o $@' Makefile
fi

echo "=== Directory Structure Created ==="
