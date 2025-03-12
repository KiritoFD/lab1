#include "../include/dna_io.h"
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

// Read DNA sequence from file - optimized for R9 7940HX
char* read_sequence_from_file(const char* filename, int* length) {
    // Use lower-level I/O for faster file reading
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        fprintf(stderr, "Could not open file: %s\n", filename);
        return NULL;
    }
    
    // Get file size more efficiently
    struct stat st;
    fstat(fd, &st);
    long file_size = st.st_size;
    
    // Use aligned memory allocation for better performance
    char* raw_buffer = allocate_dna_sequence(file_size);
    
    // Read the entire file in one operation
    ssize_t bytes_read = read(fd, raw_buffer, file_size);
    if (bytes_read != file_size) {
        fprintf(stderr, "Failed to read entire file\n");
        free(raw_buffer);
        close(fd);
        return NULL;
    }
    raw_buffer[file_size] = '\0';
    close(fd);
    
    // Pre-allocate output buffer based on file size estimate
    // Typically DNA files are mostly valid characters, so allocate at 90% of file size
    size_t estimated_seq_length = file_size * 0.9;
    char* sequence = allocate_dna_sequence(estimated_seq_length);
    
    // Process the sequence with optimized filtering
    int seq_length = 0;
    #pragma omp parallel if (file_size > 1024*1024) // Only parallelize for large files
    {
        #pragma omp for reduction(+:seq_length) schedule(static)
        for (long i = 0; i < file_size; i++) {
            char c = toupper(raw_buffer[i]);
            // Use branchless approach for character filtering
            int is_valid = (c == 'A' || c == 'T' || c == 'G' || c == 'C');
            if (is_valid) {
                // Each thread computes its local index
                int local_idx = seq_length++;
                if ((size_t)local_idx < estimated_seq_length) {
                    sequence[local_idx] = c;
                }
            }
        }
    }
    
    // Ensure proper null termination
    if ((size_t)seq_length < estimated_seq_length) {
        sequence[seq_length] = '\0';
    } else {
        // Unlikely case: resize if our estimate was too small
        sequence = realloc(sequence, seq_length + 1);
        sequence[seq_length] = '\0';
    }
    
    // Set the output length and free the raw buffer
    *length = seq_length;
    free(raw_buffer);
    
    return sequence;
}

// Save results to files
void save_results(RepeatPattern* repeats, int count, int total_repeats) {
    // Save basic results to file
    FILE* result_file = fopen("repeat_results.txt", "w");
    if (result_file) {
        fprintf(result_file, "Found Repeats:\n");
        fprintf(result_file, "Position | Length | Repeat Count | Reverse Complement\n");
        fprintf(result_file, "--------------------------------------------------\n");
        
        if (repeats) {
            for (int i = 0; i < count; i++) {
                fprintf(result_file, "%8d | %6d | %12d | %s\n", 
                       repeats[i].position, 
                       repeats[i].length, 
                       repeats[i].count, 
                       repeats[i].is_reverse ? "Yes" : "No");
            }
        }
        
        fprintf(result_file, "\nTotal: %d unique repeat positions found (total of %d repeats)\n", 
               count, total_repeats);
        
        fclose(result_file);
        printf("Basic repeat information saved to repeat_results.txt\n");
    } else {
        fprintf(stderr, "Could not open result file for writing\n");
    }
    
    // Save detailed results
    FILE* detail_file = fopen("repeat_details.txt", "w");
    if (detail_file) {
        fprintf(detail_file, "Position,Length,RepeatCount,ReverseComplement,OriginalSequence,RepeatInstances\n");
        
        if (repeats) {
            for (int i = 0; i < count; i++) {
                fprintf(detail_file, "%d,%d,%d,%s,%s,", 
                       repeats[i].position,
                       repeats[i].length,
                       repeats[i].count,
                       repeats[i].is_reverse ? "Yes" : "No",
                       repeats[i].orig_seq);
                
                // Add repeat examples
                if (repeats[i].num_examples > 0) {
                    for (int j = 0; j < repeats[i].num_examples; j++) {
                        fprintf(detail_file, "%s", repeats[i].repeat_examples[j]);
                        if (j < repeats[i].num_examples - 1) {
                            fprintf(detail_file, ";");
                        }
                    }
                } else {
                    fprintf(detail_file, "NoExamplesFound");
                }
                
                fprintf(detail_file, "\n");
            }
        }
        
        fclose(detail_file);
        printf("\nDetailed repeat information saved to repeat_details.txt\n");
    }
}
