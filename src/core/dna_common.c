#include "../include/core/dna_common.h"

// Get the reverse complement of a DNA sequence
char* get_reverse_complement(const char* sequence, int length) {
    char* result = (char*)malloc((length + 1) * sizeof(char));
    if (!result) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < length; i++) {
        char base = sequence[length - 1 - i];
        switch (base) {
            case 'A': result[i] = 'T'; break;
            case 'T': result[i] = 'A'; break;
            case 'G': result[i] = 'C'; break;
            case 'C': result[i] = 'G'; break;
            default: result[i] = 'N'; break;
        }
    }
    result[length] = '\0';
    return result;
}

// Free memory used by repeat patterns
void free_repeat_patterns(RepeatPattern* repeats, int count) {
    if (!repeats) return;
    
    for (int i = 0; i < count; i++) {
        free(repeats[i].orig_seq);
        
        if (repeats[i].repeat_examples) {
            for (int j = 0; j < repeats[i].num_examples; j++) {
                free(repeats[i].repeat_examples[j]);
            }
            free(repeats[i].repeat_examples);
        }
    }
    
    free(repeats);
}

// Print program usage instructions
void print_usage(const char* program_name) {
    printf("Usage: %s <reference_file> <query_file>\n", program_name);
    printf("Example: %s reference.txt query.txt\n", program_name);
}
