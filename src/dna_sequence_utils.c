// Utility functions for DNA sequence processing
// This replaces the old monolithic dna_repeat_finder.c file

#include "../include/dna_common.h"
#include "../include/dna_io.h"
#include "../include/dna_traditional.h"
#include "../include/dna_graph.h"

// Any additional utility functions that don't fit elsewhere can go here.
// This file serves as a placeholder replacement for dna_repeat_finder.c
// without any function definitions that would cause duplicate symbol errors.

void show_dna_stats(const char* reference, const char* query, int ref_len, int query_len) {
    int a_count = 0, t_count = 0, g_count = 0, c_count = 0;
    
    // Count nucleotides in the reference sequence
    for (int i = 0; i < ref_len; i++) {
        switch (reference[i]) {
            case 'A': a_count++; break;
            case 'T': t_count++; break;
            case 'G': g_count++; break;
            case 'C': c_count++; break;
            default: break; // Should not happen with filtered input
        }
    }
    
    printf("\nReference Sequence Statistics:\n");
    printf("Length: %d bases\n", ref_len);
    printf("A: %d (%.1f%%)\n", a_count, (float)a_count * 100 / ref_len);
    printf("T: %d (%.1f%%)\n", t_count, (float)t_count * 100 / ref_len);
    printf("G: %d (%.1f%%)\n", g_count, (float)g_count * 100 / ref_len);
    printf("C: %d (%.1f%%)\n", c_count, (float)c_count * 100 / ref_len);
    
    // Reset counters for query sequence
    a_count = 0;
    t_count = 0;
    g_count = 0;
    c_count = 0;
    
    // Count nucleotides in the query sequence
    for (int i = 0; i < query_len; i++) {
        switch (query[i]) {
            case 'A': a_count++; break;
            case 'T': t_count++; break;
            case 'G': g_count++; break;
            case 'C': c_count++; break;
            default: break; // Should not happen with filtered input
        }
    }
    
    printf("\nQuery Sequence Statistics:\n");
    printf("Length: %d bases\n", query_len);
    printf("A: %d (%.1f%%)\n", a_count, (float)a_count * 100 / query_len);
    printf("T: %d (%.1f%%)\n", t_count, (float)t_count * 100 / query_len);
    printf("G: %d (%.1f%%)\n", g_count, (float)g_count * 100 / query_len);
    printf("C: %d (%.1f%%)\n", c_count, (float)c_count * 100 / query_len);
}
