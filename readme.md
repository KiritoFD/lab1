# DNA Repeat Finder

A high-performance application for finding repeat patterns in DNA sequences, optimized for AMD Ryzen 9 7940HX processors.

## Project Overview

This application implements two approaches for finding DNA repeats:
1. **Traditional Method**: Direct string matching with parallel optimization
2. **DAG-based Method**: Models DNA sequences as a directed acyclic graph for path planning

## Project Structure

```
/home/xy/lab1/
|-- src/                  # Source code files
|   |-- main.c            # Main program entry point
|   |-- dna_common.c      # Common DNA utilities
|   |-- dna_io.c          # I/O operations for DNA sequences
|   |-- dna_traditional.c # Traditional repeat finding implementation
|   |-- dna_graph.c       # Graph-based repeat finding implementation
|
|-- include/              # Header files
|   |-- dna_common.h      # Common DNA utilities declarations
|   |-- dna_io.h          # I/O operations declarations
|   |-- dna_traditional.h # Traditional approach declarations
|   |-- dna_graph.h       # Graph-based approach declarations
|   |-- cpu_optimize.h    # R9 7940HX specific optimizations
|
|-- bin/                  # Executable files (created during build)
|-- obj/                  # Object files (created during build)
|-- Makefile              # Build configuration
|-- README.md             # This file
```

## Requirements

- GCC compiler with C11 support
- OpenMP for parallel processing
- Make build system

## Compilation

The project includes a Makefile with various build targets optimized for the R9 7940HX processor.

### Basic build:
```bash
cd /home/xy/lab1
make
```

This will automatically use all available CPU cores for compilation and apply optimizations for the R9 7940HX.

### Other build options:

```bash
# Build with debugging information
make debug

# Clean compiled objects
make clean

# Complete clean (including binaries)
make distclean

# Display help
make help
```

## Running the Program

### Basic usage:
```bash
./bin/dna_repeat_finder reference_file query_file
```

Example:
```bash
./bin/dna_repeat_finder sample_data/reference.txt sample_data/query.txt
```

### Input Files

The program accepts two text files:
1. **Reference file**: Contains the reference DNA sequence
2. **Query file**: Contains the DNA sequence to search for repeats

Both files should contain DNA sequences (A, T, G, C). Non-DNA characters are automatically filtered out.

### Output Files

Two output files are generated:
1. **repeat_results.txt**: Basic information about found repeats
2. **repeat_details.txt**: Detailed information including sequence examples

## Performance Optimizations

This version includes extensive optimizations for the AMD Ryzen 9 7940HX processor:

- Zen 4 architecture-specific compiler flags
- AVX2 vectorization for DNA comparison operations
- Cache-aligned memory allocation
- Memory prefetching for optimal cache utilization
- Dynamic thread allocation based on workload size
- Link-time optimization (LTO)
- Loop optimization with auto-vectorization

## Algorithm Description

### Traditional Approach
The traditional method searches for repeats by:
1. Extracting segments from the reference sequence
2. Searching for these segments in the query sequence
3. Checking for consecutive repeats
4. Also checking reverse complement matches

### Graph-based DAG Approach
The DAG approach models the problem as:
1. Nodes represent positions in the DNA sequence
2. Edges represent matching segments between positions
3. Path planning algorithms find optimal paths through the graph
4. These paths correspond to repeat patterns

## License

[Your License Information]

## Authors

[Your Author Information]