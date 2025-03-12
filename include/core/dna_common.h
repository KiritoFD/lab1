#ifndef DNA_COMMON_H
#define DNA_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <omp.h>
#include "cpu_optimize.h" // Include R9 7940HX specific optimizations

// Structure to represent a repeat pattern
typedef struct {
    int position;       // Position in reference sequence
    int length;         // Length of the repeat
    int count;          // Number of repeats
    int is_reverse;     // 1 if reverse complement, 0 otherwise
    char* orig_seq ALIGN_TO_CACHE; // Original sequence, cache-aligned
    char** repeat_examples; // Examples of repeat sequences found
    int num_examples;   // Number of examples stored
} RepeatPattern;

// Common utility functions
char* get_reverse_complement(const char* sequence, int length);
void free_repeat_patterns(RepeatPattern* repeats, int count);
void print_usage(const char* program_name);

// Optimized memory allocation for DNA sequences
FORCE_INLINE char* allocate_dna_sequence(size_t length) {
    // Use aligned allocation for better memory access performance
    char* sequence = (char*)aligned_alloc(CACHE_LINE_SIZE, ((length + 1 + CACHE_LINE_SIZE-1) & ~(CACHE_LINE_SIZE-1)));
    if (UNLIKELY(!sequence)) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    return sequence;
}

#endif // DNA_COMMON_H
