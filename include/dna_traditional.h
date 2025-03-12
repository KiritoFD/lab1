#ifndef DNA_TRADITIONAL_H
#define DNA_TRADITIONAL_H

#include "dna_common.h"

// Functions for traditional approach
int** build_similarity_matrix(const char* reference, int ref_len, const char* query, int query_len);
RepeatPattern* find_repeats(const char* reference, int ref_len, const char* query, int query_len, int* num_repeats);
RepeatPattern* get_repeat_sequences(RepeatPattern* repeats, int num_repeats, const char* reference, const char* query, int ref_len, int query_len, int* new_count);
RepeatPattern* filter_nested_repeats(RepeatPattern* repeats, int num_repeats, int filter_no_instances, int* filtered_count);
void free_matrix(int** matrix, int rows);

#endif // DNA_TRADITIONAL_H
