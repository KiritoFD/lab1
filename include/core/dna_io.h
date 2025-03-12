#ifndef DNA_IO_H
#define DNA_IO_H

#include "dna_common.h"

// File and sequence handling functions
char* read_sequence_from_file(const char* filename, int* length);
void save_results(RepeatPattern* repeats, int count, int total_repeats);

#endif // DNA_IO_H
