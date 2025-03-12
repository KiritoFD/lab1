#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

// Generate a random DNA sequence of specified length
char* generate_dna_sequence(int length) {
    char* sequence = (char*)malloc((length + 1) * sizeof(char));
    if (!sequence) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    const char bases[] = "ACGT";
    
    // Generate sequence in chunks to improve performance
    int chunk_size = 10000;
    for (int i = 0; i < length; i += chunk_size) {
        int current_chunk = (i + chunk_size < length) ? chunk_size : (length - i);
        for (int j = 0; j < current_chunk; j++) {
            sequence[i + j] = bases[rand() % 4];
        }
    }
    sequence[length] = '\0';
    
    return sequence;
}

// Insert repeats into a DNA sequence
void insert_repeats(char* sequence, int seq_length) {
    // Number of repeats scales with sequence length
    int num_repeats = seq_length / 10000;  // One repeat per 10k bases
    
    printf("Inserting %d repeat patterns...\n", num_repeats);
    
    // Short repeats (5-15 bp)
    for (int i = 0; i < num_repeats; i++) {
        int repeat_len = 5 + rand() % 11;  // 5 to 15 bases
        char* repeat = (char*)malloc((repeat_len + 1) * sizeof(char));
        
        // Generate random repeat
        for (int j = 0; j < repeat_len; j++) {
            repeat[j] = "ACGT"[rand() % 4];
        }
        repeat[repeat_len] = '\0';
        
        // Insert at 2-4 locations
        int num_locations = 2 + rand() % 3;
        for (int j = 0; j < num_locations; j++) {
            int pos = rand() % (seq_length - repeat_len - 1);
            strncpy(sequence + pos, repeat, repeat_len);
        }
        
        free(repeat);
    }
    
    // Medium repeats (20-50 bp)
    for (int i = 0; i < num_repeats / 2; i++) {
        int repeat_len = 20 + rand() % 31;  // 20 to 50 bases
        char* repeat = (char*)malloc((repeat_len + 1) * sizeof(char));
        
        // Generate random repeat
        for (int j = 0; j < repeat_len; j++) {
            repeat[j] = "ACGT"[rand() % 4];
        }
        repeat[repeat_len] = '\0';
        
        // Insert at 1-3 locations
        int num_locations = 1 + rand() % 3;
        for (int j = 0; j < num_locations; j++) {
            int pos = rand() % (seq_length - repeat_len - 1);
            strncpy(sequence + pos, repeat, repeat_len);
        }
        
        free(repeat);
    }
    
    // Long repeats (80-120 bp)
    for (int i = 0; i < num_repeats / 5; i++) {
        int repeat_len = 80 + rand() % 41;  // 80 to 120 bases
        char* repeat = (char*)malloc((repeat_len + 1) * sizeof(char));
        
        // Generate random repeat
        for (int j = 0; j < repeat_len; j++) {
            repeat[j] = "ACGT"[rand() % 4];
        }
        repeat[repeat_len] = '\0';
        
        // Insert at 1-2 locations
        int num_locations = 1 + rand() % 2;
        for (int j = 0; j < num_locations; j++) {
            int pos = rand() % (seq_length - repeat_len - 1);
            strncpy(sequence + pos, repeat, repeat_len);
        }
        
        free(repeat);
    }
    
    // Insert some reverse complement repeats
    for (int i = 0; i < num_repeats / 2; i++) {
        int repeat_len = 10 + rand() % 41;  // 10 to 50 bases
        char* repeat = (char*)malloc((repeat_len + 1) * sizeof(char));
        char* rev_comp = (char*)malloc((repeat_len + 1) * sizeof(char));
        
        // Generate random repeat
        for (int j = 0; j < repeat_len; j++) {
            repeat[j] = "ACGT"[rand() % 4];
        }
        repeat[repeat_len] = '\0';
        
        // Create reverse complement
        for (int j = 0; j < repeat_len; j++) {
            char base = repeat[repeat_len - 1 - j];
            switch (base) {
                case 'A': rev_comp[j] = 'T'; break;
                case 'T': rev_comp[j] = 'A'; break;
                case 'G': rev_comp[j] = 'C'; break;
                case 'C': rev_comp[j] = 'G'; break;
                default: rev_comp[j] = 'N'; break;
            }
        }
        rev_comp[repeat_len] = '\0';
        
        // Insert original and reverse complement
        int pos1 = rand() % (seq_length / 2 - repeat_len - 1);
        int pos2 = (seq_length / 2) + rand() % (seq_length / 2 - repeat_len - 1);
        
        strncpy(sequence + pos1, repeat, repeat_len);
        strncpy(sequence + pos2, rev_comp, repeat_len);
        
        free(repeat);
        free(rev_comp);
    }
}

// Write sequence to file efficiently
void write_sequence_to_file(FILE* file, const char* sequence, int length, const char* header) {
    fprintf(file, ">%s\n", header);
    
    // Write in chunks of 80 characters per line
    const int line_length = 80;
    for (int i = 0; i < length; i += line_length) {
        int current_line = (i + line_length < length) ? line_length : (length - i);
        fprintf(file, "%.*s\n", current_line, sequence + i);
    }
}

int main(int argc, char* argv[]) {
    // Seed the random number generator
    srand((unsigned int)time(NULL));
    
    // Default settings
    int length = 100000;
    char ref_filename[256] = "reference.txt";
    char query_filename[256] = "query.txt";
    
    // Process command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-length") == 0 && i + 1 < argc) {
            length = atoi(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "-ref") == 0 && i + 1 < argc) {
            strncpy(ref_filename, argv[i + 1], sizeof(ref_filename) - 1);
            i++;
        } else if (strcmp(argv[i], "-query") == 0 && i + 1 < argc) {
            strncpy(query_filename, argv[i + 1], sizeof(query_filename) - 1);
            i++;
        }
    }
    
    printf("Generating DNA sequences of length %d...\n", length);
    
    // Generate reference sequence
    clock_t start_time = clock();
    printf("Generating reference sequence...\n");
    char* reference = generate_dna_sequence(length);
    printf("Reference sequence generated\n");
    
    // Generate query sequence (a modified copy of reference)
    printf("Generating query sequence...\n");
    char* query = (char*)malloc((length + 1) * sizeof(char));
    if (!query) {
        fprintf(stderr, "Memory allocation failed\n");
        free(reference);
        exit(EXIT_FAILURE);
    }
    
    // Copy reference to query
    memcpy(query, reference, length + 1);
    printf("Query sequence created\n");
    
    // Insert repeats in both sequences
    printf("Inserting repeats in reference sequence...\n");
    insert_repeats(reference, length);
    printf("Inserting repeats in query sequence...\n");
    insert_repeats(query, length);
    
    // Add some differences between reference and query (about 5% difference)
    printf("Adding mutations between sequences...\n");
    int num_differences = length / 20;
    for (int i = 0; i < num_differences; i++) {
        int pos = rand() % length;
        char new_base;
        do {
            new_base = "ACGT"[rand() % 4];
        } while (new_base == query[pos]);
        query[pos] = new_base;
    }
    
    // Calculate time taken
    clock_t end_time = clock();
    double generation_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Sequence generation completed in %.2f seconds\n", generation_time);
    
    // Save to files
    printf("Writing reference sequence to %s...\n", ref_filename);
    FILE* ref_file = fopen(ref_filename, "w");
    if (!ref_file) {
        fprintf(stderr, "Could not open file for writing: %s\n", ref_filename);
        free(reference);
        free(query);
        exit(EXIT_FAILURE);
    }
    
    write_sequence_to_file(ref_file, reference, length, "Reference sequence length 100000");
    fclose(ref_file);
    
    printf("Writing query sequence to %s...\n", query_filename);
    FILE* query_file = fopen(query_filename, "w");
    if (!query_file) {
        fprintf(stderr, "Could not open file for writing: %s\n", query_filename);
        free(reference);
        free(query);
        exit(EXIT_FAILURE);
    }
    
    write_sequence_to_file(query_file, query, length, "Query sequence length 100000");
    fclose(query_file);
    
    printf("\nGenerated sequences successfully:\n");
    printf("- Reference sequence saved to: %s\n", ref_filename);
    printf("- Query sequence saved to: %s\n", query_filename);
    printf("- Length: %d base pairs\n", length);
    printf("- Generation time: %.2f seconds\n", generation_time);
    
    // Clean up
    free(reference);
    free(query);
    
    return 0;
}
