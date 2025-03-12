# DNA Repeat Finder

A high-performance application for finding repeat patterns in DNA sequences, optimized for AMD Ryzen 9 7940HX processors.

## Project Structure

```
/home/xy/lab1/
├── src/               # Core C implementation
│   ├── main.c         # Main program entry point
│   ├── dna_common.c   # Common DNA utilities
│   ├── dna_io.c       # I/O operations for DNA sequences
│   ├── dna_traditional.c # Traditional repeat finding implementation
│   ├── dna_graph.c    # Graph-based repeat finding implementation
│   └── dna_sequence_utils.c # Sequence utilities
│
├── include/           # Header files
│   ├── dna_common.h   # Common DNA utilities declarations
│   ├── dna_io.h       # I/O operations declarations
│   ├── dna_traditional.h # Traditional approach declarations
│   ├── dna_graph.h    # Graph-based approach declarations
│   └── cpu_optimize.h # R9 7940HX specific optimizations
│
├── tools/            # Additional tools
│   └── generate_test_sequences.c # Test sequence generator (optional)
│
├── algorithms/       # DAG algorithm implementations (optional)
│   ├── dag_algorithms.c
│   ├── dag_algorithms.h
│   └── dag_path_planning.md
│
├── cpp/             # C++ implementation (alternative)
│   └── dna_repeat_finder.cpp
│
├── bin/             # Executable files (created during build)
├── obj/             # Object files (created during build)
├── Makefile         # Build configuration
└── README.md        # This file
```

## Building the Project

```bash
# Standard build
make

# Debug build
make debug

# Clean and rebuild
make clean && make
```

## Running the Program

```bash
# Using default files (reference.txt and query.txt)
./bin/dna_repeat_finder

# Specifying input files
./bin/dna_repeat_finder path/to/reference.txt path/to/query.txt
```
