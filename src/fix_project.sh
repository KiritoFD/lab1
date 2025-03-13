#!/bin/bash

echo "=== DNA Repeat Finder Project Structure Fix ==="

# Create necessary directories if they don't exist
echo "Creating project directories..."
mkdir -p bin obj data

# Clean up object files
echo "Cleaning object files..."
rm -rf obj/*.o obj/*.d

# Fix duplicate main function issue
if [ -f "src/generate_test_sequences.c" ]; then
    echo "Fixing duplicate main function in generate_test_sequences.c..."
    cp src/generate_test_sequences.c src/generate_test_sequences.c.backup
    sed -i 's/int main(int argc, char\* argv\[\]) {/int generator_main(int argc, char* argv[]) {/' src/generate_test_sequences.c
    
    # Add conditional main function at the end of the file
    cat >> src/generate_test_sequences.c << 'EOF'

#ifdef GENERATOR_BUILD
// Only include this main function when specifically building the generator
int main(int argc, char* argv[]) {
    return generator_main(argc, argv);
}
#endif
EOF
    echo "  - Modified generate_test_sequences.c to use generator_main"
fi

# Fix Makefile to handle conditional compilation
echo "Updating Makefile..."
if [ -f "Makefile" ]; then
    cp Makefile Makefile.backup
    
    # Add conditional compilation for the generator
    sed -i 's|$(GENERATOR_TARGET): $(GENERATOR_OBJECT)|$(GENERATOR_TARGET): $(COMMON_OBJECTS) $(GENERATOR_OBJECT)|' Makefile
    sed -i 's|$(CC) $^ -o $@ $(LDFLAGS)|$(CC) $(CFLAGS) -DGENERATOR_BUILD $^ -o $@ $(LDFLAGS)|' Makefile
    
    echo "  - Updated Makefile for conditional compilation"
fi

# Fix include paths if needed
echo "Checking include paths..."
if grep -q "../include/core/" src/main.c; then
    echo "  - Fixing incorrect include paths in main.c"
    sed -i 's|../include/core/|../include/|g' src/main.c
fi

# Make sure src/dna_graph.h includes are correctly pointing to the main include directory
if [ -f "src/dna_graph.h" ]; then
    echo "Moving src/dna_graph.h to include directory if needed..."
    if ! [ -f "include/dna_graph.h" ]; then
        cp src/dna_graph.h include/
        echo "  - Copied src/dna_graph.h to include/"
    fi
fi

# Check for dna_common.h location and fix
if [ -f "src/dna_common.h" ] && ! [ -f "include/dna_common.h" ]; then
    echo "Fixing dna_common.h location..."
    cp src/dna_common.h include/
    echo "  - Copied src/dna_common.h to include/"
fi

# Fix dna_io.h location
if [ -f "src/dna_io.h" ] && ! [ -f "include/dna_io.h" ]; then
    echo "Fixing dna_io.h location..."
    cp src/dna_io.h include/
    echo "  - Copied src/dna_io.h to include/"
fi

# Fix permissions
echo "Setting executable permissions..."
chmod +x fix_project.sh

# Create default reference and query files if they don't exist
if [ ! -f "reference.txt" ]; then
    echo "Creating default reference.txt..."
    touch reference.txt
fi

if [ ! -f "query.txt" ]; then
    echo "Creating default query.txt..."
    touch query.txt
fi

echo ""
echo "Project structure fixed. Compile with:"
echo "  make clean && make"
echo ""
echo "=== Done ==="
