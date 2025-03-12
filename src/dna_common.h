#ifndef DNA_COMMON_H
#define DNA_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <omp.h>

// Structure to represent a repeat pattern
typedef struct {
    int position;       // Position in reference sequence
    int length;         // Length of the repeat
    int count;          // Number of repeats
    int is_reverse;     // 1 if reverse complement, 0 otherwise
    char* orig_seq;     // Original sequence
    char** repeat_examples; // Examples of repeat sequences found
    int num_examples;   // Number of examples stored
} RepeatPattern;

// Common utility functions
char* get_reverse_complement(const char* sequence, int length);
void free_repeat_patterns(RepeatPattern* repeats, int count);
void print_usage(const char* program_name);

#endif // DNA_COMMON_H
