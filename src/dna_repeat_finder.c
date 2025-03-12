#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <omp.h> // Add OpenMP support

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

// Function prototypes
char* get_reverse_complement(const char* sequence, int length);
int** build_similarity_matrix(const char* reference, int ref_len, const char* query, int query_len);
RepeatPattern* find_repeats(const char* reference, int ref_len, const char* query, int query_len, int* num_repeats);
RepeatPattern* get_repeat_sequences(RepeatPattern* repeats, int num_repeats, const char* reference, const char* query, int ref_len, int query_len, int* new_count);
RepeatPattern* filter_nested_repeats(RepeatPattern* repeats, int num_repeats, int filter_no_instances, int* filtered_count);
void free_repeat_patterns(RepeatPattern* repeats, int count);
void free_matrix(int** matrix, int rows);

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

// Build similarity matrix between reference and query - optimized with parallel processing
int** build_similarity_matrix(const char* reference, int ref_len, const char* query, int query_len) {
    int** matrix = (int**)malloc(ref_len * sizeof(int*));
    if (!matrix) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < ref_len; i++) {
        matrix[i] = (int*)malloc(query_len * sizeof(int));
        if (!matrix[i]) {
            fprintf(stderr, "Memory allocation failed in thread %d\n", omp_get_thread_num());
            // Cannot safely free memory here due to parallel context
            exit(EXIT_FAILURE);
        }
        
        for (int j = 0; j < query_len; j++) {
            matrix[i][j] = (reference[i] == query[j]) ? 1 : -1;
        }
    }
    
    return matrix;
}

// Find repeats in the query sequence compared to reference
RepeatPattern* find_repeats(const char* reference, int ref_len, const char* query, int query_len, int* num_repeats) {
    // Dynamic parameters based on sequence length
    int min_length = 5 > (ref_len / 1000) ? 5 : (ref_len / 1000);
    int max_length = (ref_len / 10) < 120 ? (ref_len / 10) : 120;
    int step = 1 > (min_length / 5) ? 1 : (min_length / 5);
    
    printf("Using parameters: min_length=%d, max_length=%d, step=%d\n", 
           min_length, max_length, step);
    
    // Limit search space for very large sequences
    int max_positions_to_check = 10000;
    int positions_step = ref_len > max_positions_to_check ? (ref_len / max_positions_to_check) : 1;
    
    // Allocate initial memory for repeats (will be reallocated as needed)
    int capacity = 1000;
    RepeatPattern* repeats = (RepeatPattern*)malloc(capacity * sizeof(RepeatPattern));
    if (!repeats) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    int repeat_count = 0;
    
    // Progress reporting
    int progress = 0;
    int total_progress = (ref_len / positions_step);
    int last_reported_progress = -1;
    
    // Check positions in the reference sequence with step size for large sequences
    for (int pos = 0; pos < ref_len - min_length; pos += positions_step) {
        // Report progress
        progress++;
        int current_progress = (progress * 100) / total_progress;
        if (current_progress % 10 == 0 && current_progress != last_reported_progress) {
            printf("Finding repeats: %d%% complete\n", current_progress);
            last_reported_progress = current_progress;
        }
        
        // Check different segment lengths
        for (int length = min_length; length < (max_length < (ref_len - pos) ? max_length : (ref_len - pos)); length += step) {
            // Extract segment from reference
            char* segment = (char*)malloc((length + 1) * sizeof(char));
            if (!segment) {
                fprintf(stderr, "Memory allocation failed\n");
                free_repeat_patterns(repeats, repeat_count);
                exit(EXIT_FAILURE);
            }
            
            strncpy(segment, reference + pos, length);
            segment[length] = '\0';
            
            // Check for forward repeats
            int start_idx = 0;
            while (1) {
                // Find next occurrence of segment
                char* found = strstr(query + start_idx, segment);
                if (!found) break;
                
                int next_idx = found - query;
                int current_pos = next_idx + length;
                int consecutive_count = 0;
                
                // Check for consecutive repeats
                while (current_pos + length <= query_len) {
                    int is_repeat = 1;
                    for (int k = 0; k < length; k++) {
                        if (query[current_pos + k] != segment[k]) {
                            is_repeat = 0;
                            break;
                        }
                    }
                    
                    if (is_repeat) {
                        consecutive_count++;
                        current_pos += length;
                    } else {
                        break;
                    }
                }
                
                if (consecutive_count > 0) {
                    // Add to repeats array
                    if (repeat_count >= capacity) {
                        capacity *= 2;
                        RepeatPattern* new_repeats = (RepeatPattern*)realloc(repeats, capacity * sizeof(RepeatPattern));
                        if (!new_repeats) {
                            fprintf(stderr, "Memory reallocation failed\n");
                            free(segment);
                            free_repeat_patterns(repeats, repeat_count);
                            exit(EXIT_FAILURE);
                        }
                        repeats = new_repeats;
                    }
                    
                    repeats[repeat_count].position = pos;
                    repeats[repeat_count].length = length;
                    repeats[repeat_count].count = consecutive_count;
                    repeats[repeat_count].is_reverse = 0;
                    repeats[repeat_count].orig_seq = strdup(segment);
                    repeats[repeat_count].repeat_examples = NULL;
                    repeats[repeat_count].num_examples = 0;
                    repeat_count++;
                }
                
                start_idx = next_idx + 1;
                if (start_idx >= query_len) break;
            }
            
            // Check for reverse complement repeats
            char* rev_comp = get_reverse_complement(segment, length);
            start_idx = 0;
            
            while (1) {
                char* found = strstr(query + start_idx, rev_comp);
                if (!found) break;
                
                int next_idx = found - query;
                int current_pos = next_idx + length;
                int consecutive_count = 0;
                
                // Check for consecutive reverse repeats
                while (current_pos + length <= query_len) {
                    int is_repeat = 1;
                    for (int k = 0; k < length; k++) {
                        if (query[current_pos + k] != rev_comp[k]) {
                            is_repeat = 0;
                            break;
                        }
                    }
                    
                    if (is_repeat) {
                        consecutive_count++;
                        current_pos += length;
                    } else {
                        break;
                    }
                }
                
                // Add to repeats array (always add for reverse complement)
                if (repeat_count >= capacity) {
                    capacity *= 2;
                    RepeatPattern* new_repeats = (RepeatPattern*)realloc(repeats, capacity * sizeof(RepeatPattern));
                    if (!new_repeats) {
                        fprintf(stderr, "Memory reallocation failed\n");
                        free(segment);
                        free(rev_comp);
                        free_repeat_patterns(repeats, repeat_count);
                        exit(EXIT_FAILURE);
                    }
                    repeats = new_repeats;
                }
                
                repeats[repeat_count].position = pos;
                repeats[repeat_count].length = length;
                repeats[repeat_count].count = consecutive_count > 0 ? consecutive_count : 1;
                repeats[repeat_count].is_reverse = 1;
                repeats[repeat_count].orig_seq = strdup(segment);
                repeats[repeat_count].repeat_examples = NULL;
                repeats[repeat_count].num_examples = 0;
                repeat_count++;
                
                start_idx = next_idx + 1;
                if (start_idx >= query_len) break;
            }
            
            free(segment);
            free(rev_comp);
        }
    }
    
    printf("Found %d repeat patterns\n", repeat_count);
    *num_repeats = repeat_count;
    
    // If no repeats found, free memory and return NULL
    if (repeat_count == 0) {
        free(repeats);
        return NULL;
    }
    
    return repeats;
}

// Add sequence information to repeat patterns - parallelized version
RepeatPattern* get_repeat_sequences(RepeatPattern* repeats, int num_repeats, const char* reference, const char* query, int ref_len, int query_len, int* new_count) {
    *new_count = num_repeats;
    
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num_repeats; i++) {
        int pos = repeats[i].position;
        int length = repeats[i].length;
        int is_reverse = repeats[i].is_reverse;
        
        // Get the original segment from reference
        char* segment = (char*)malloc((length + 1) * sizeof(char));
        if (!segment) {
            fprintf(stderr, "Memory allocation failed in thread %d\n", omp_get_thread_num());
            exit(EXIT_FAILURE);
        }
        strncpy(segment, reference + pos, length);
        segment[length] = '\0';
        
        // For reverse complement
        char* search_seq;
        if (is_reverse) {
            search_seq = get_reverse_complement(segment, length);
        } else {
            search_seq = strdup(segment);
        }
        
        // Find repeat instances in query
        int start_idx = 0;
        int instances_count = 0;
        int instances_capacity = 10;
        char** instances = (char**)malloc(instances_capacity * sizeof(char*));
        
        while (1) {
            char* found = strstr(query + start_idx, search_seq);
            if (!found) break;
            
            int next_idx = found - query;
            
            if (instances_count >= instances_capacity) {
                instances_capacity *= 2;
                char** new_instances = (char**)realloc(instances, instances_capacity * sizeof(char*));
                if (!new_instances) {
                    fprintf(stderr, "Memory reallocation failed in thread %d\n", omp_get_thread_num());
                    free(segment);
                    free(search_seq);
                    for (int j = 0; j < instances_count; j++) {
                        free(instances[j]);
                    }
                    free(instances);
                    exit(EXIT_FAILURE);
                }
                instances = new_instances;
            }
            
            instances[instances_count] = (char*)malloc((length + 1) * sizeof(char));
            strncpy(instances[instances_count], query + next_idx, length);
            instances[instances_count][length] = '\0';
            instances_count++;
            
            start_idx = next_idx + 1;
            if (start_idx >= query_len) break;
            if (instances_count >= repeats[i].count) break; // Only collect up to count examples
        }
        
        #pragma omp critical
        {
            repeats[i].repeat_examples = instances;
            repeats[i].num_examples = instances_count;
        }
        
        // Free temporary memory
        free(search_seq);
    }
    
    return repeats;
}

// Filter nested repeats, keeping only the longest at each position
RepeatPattern* filter_nested_repeats(RepeatPattern* repeats, int num_repeats, int filter_no_instances, int* filtered_count) {
    if (num_repeats == 0) {
        *filtered_count = 0;
        return NULL;
    }
    
    // First filter by instances if requested
    RepeatPattern* filtered = repeats;
    int filtered_size = num_repeats;
    
    if (filter_no_instances) {
        int count = 0;
        RepeatPattern* temp = (RepeatPattern*)malloc(num_repeats * sizeof(RepeatPattern));
        
        for (int i = 0; i < num_repeats; i++) {
            if (repeats[i].num_examples > 0) {
                temp[count++] = repeats[i];
            }
        }
        
        // If we filtered out everything, revert to original
        if (count == 0) {
            free(temp);
            filtered = repeats;
            filtered_size = num_repeats;
        } else {
            filtered = temp;
            filtered_size = count;
        }
    }
    
    // Create a temporary array to track positions already processed
    typedef struct {
        int position;
        int is_reverse;
        int used;
    } PositionKey;
    
    PositionKey* positions = (PositionKey*)malloc(filtered_size * sizeof(PositionKey));
    
    for (int i = 0; i < filtered_size; i++) {
        positions[i].position = filtered[i].position;
        positions[i].is_reverse = filtered[i].is_reverse;
        positions[i].used = 0;
    }
    
    // Allocate result array
    RepeatPattern* result = (RepeatPattern*)malloc(filtered_size * sizeof(RepeatPattern));
    int result_count = 0;
    
    // For each unique position, find the longest repeat
    for (int i = 0; i < filtered_size; i++) {
        if (positions[i].used) continue;
        
        int best_idx = i;
        int best_length = filtered[i].length;
        
        // Find the longest repeat at this position
        for (int j = i + 1; j < filtered_size; j++) {
            if (positions[j].position == positions[i].position && 
                positions[j].is_reverse == positions[i].is_reverse) {
                positions[j].used = 1; // Mark as processed
                
                if (filtered[j].length > best_length) {
                    best_idx = j;
                    best_length = filtered[j].length;
                }
            }
        }
        
        // Copy the best repeat to result
        result[result_count++] = filtered[best_idx];
    }
    
    free(positions);
    
    // If we created a temporary filtered array, free it
    if (filtered != repeats) {
        free(filtered);
    }
    
    *filtered_count = result_count;
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

// Free memory used by similarity matrix
void free_matrix(int** matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Read DNA sequence from file
char* read_sequence_from_file(const char* filename, int* length) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Could not open file: %s\n", filename);
        return NULL;
    }
    
    // Get file size
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);
    
    // Allocate memory for sequence (plus space for null terminator)
    char* sequence = (char*)malloc((file_size + 1) * sizeof(char));
    if (!sequence) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }
    
    // Read sequence and filter out non-DNA characters
    int seq_length = 0;
    int c;
    while ((c = fgetc(file)) != EOF) {
        c = toupper(c);
        if (c == 'A' || c == 'T' || c == 'G' || c == 'C') {
            sequence[seq_length++] = c;
        }
    }
    
    sequence[seq_length] = '\0';
    *length = seq_length;
    
    fclose(file);
    return sequence;
}

void print_usage(const char* program_name) {
    printf("Usage: %s <reference_file> <query_file>\n", program_name);
    printf("Example: %s reference.txt query.txt\n", program_name);
}

// Main function
int main(int argc, char *argv[]) {
    // Process command line arguments
    if (argc != 3) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    
    char* reference_file = argv[1];
    char* query_file = argv[2];
    
    printf("DNA Repeat Finder\n");
    printf("Version: 1.1 with improved error handling\n\n");
    
    // Thread handling (if OpenMP is enabled)
    #ifdef _OPENMP
    int max_threads = omp_get_max_threads();
    int thread_count = max_threads; // Default to maximum
    
    // Automatically use maximum threads without asking
    printf("Using %d threads for parallel processing\n", thread_count);
    omp_set_num_threads(thread_count);
    #endif
    
    char* reference = NULL;
    char* query = NULL;
    int ref_len = 0;
    int query_len = 0;
    
    // Read sequences from the provided files
    printf("Reading reference sequence from %s...\n", reference_file);
    reference = read_sequence_from_file(reference_file, &ref_len);
    if (!reference) {
        fprintf(stderr, "Failed to read reference file: %s\n", reference_file);
        return EXIT_FAILURE;
    }
    
    printf("Reading query sequence from %s...\n", query_file);
    query = read_sequence_from_file(query_file, &query_len);
    if (!query) {
        fprintf(stderr, "Failed to read query file: %s\n", query_file);
        free(reference);
        return EXIT_FAILURE;
    }
    
    printf("Successfully loaded sequences. Reference length: %d, Query length: %d\n", ref_len, query_len);
    
    // Record start time
    clock_t start_time = clock();
    
    // Find repeats with error handling
    int num_repeats = 0;
    RepeatPattern* repeats = NULL;
    
    try {
        repeats = find_repeats(reference, ref_len, query, query_len, &num_repeats);
        if (!repeats && num_repeats > 0) {
            fprintf(stderr, "Error: Failed to allocate memory for repeats\n");
            free(reference);
            free(query);
            return EXIT_FAILURE;
        }
    } catch (...) {
        fprintf(stderr, "Error: Exception during repeat finding\n");
        free(reference);
        free(query);
        return EXIT_FAILURE;
    }
    
    // Record end time
    clock_t end_time = clock();
    
    // Add sequence information to repeats
    int num_with_sequences = 0;
    RepeatPattern* repeats_with_sequences = NULL;
    
    if (repeats) {
        repeats_with_sequences = get_repeat_sequences(repeats, num_repeats, 
                                                     reference, query, 
                                                     ref_len, query_len,
                                                     &num_with_sequences);
    }
    
    // Filter nested repeats
    int filtered_count = 0;
    RepeatPattern* filtered_repeats = NULL;
    
    if (repeats_with_sequences) {
        filtered_repeats = filter_nested_repeats(repeats_with_sequences, num_with_sequences, 
                                               1, &filtered_count);
    }
    
    // Save basic results to file
    FILE* result_file = fopen("repeat_results.txt", "w");
    if (result_file) {
        fprintf(result_file, "Found Repeats:\n");
        fprintf(result_file, "Position | Length | Repeat Count | Reverse Complement\n");
        fprintf(result_file, "--------------------------------------------------\n");
        
        if (filtered_repeats) {
            for (int i = 0; i < filtered_count; i++) {
                fprintf(result_file, "%8d | %6d | %12d | %s\n", 
                       filtered_repeats[i].position, 
                       filtered_repeats[i].length, 
                       filtered_repeats[i].count, 
                       filtered_repeats[i].is_reverse ? "Yes" : "No");
            }
        }
        
        fprintf(result_file, "\nTotal: %d unique repeat positions found (total of %d repeats)\n", 
               filtered_count, num_repeats);
        
        fclose(result_file);
        printf("Basic repeat information saved to repeat_results.txt\n");
    } else {
        fprintf(stderr, "Could not open result file for writing\n");
    }
    
    // Calculate execution time
    double execution_time = ((double)(end_time - start_time) * 1000.0) / CLOCKS_PER_SEC;
    printf("Repeat search time: %.2f milliseconds\n", execution_time);
    
    // Display results
    printf("\nFound Repeats:\n");
    printf("Position | Length | Repeat Count | Reverse Complement\n");
    printf("--------------------------------------------------\n");
    if (filtered_repeats) {
        for (int i = 0; i < filtered_count; i++) {
            printf("%8d | %6d | %12d | %s\n", 
                  filtered_repeats[i].position, 
                  filtered_repeats[i].length, 
                  filtered_repeats[i].count, 
                  filtered_repeats[i].is_reverse ? "Yes" : "No");
        }
    }
    
    printf("\nTotal: %d unique repeat positions found (total of %d repeats)\n", filtered_count, num_repeats);
    
    // Save detailed results
    FILE* detail_file = fopen("repeat_details.txt", "w");
    if (detail_file) {
        fprintf(detail_file, "Position,Length,RepeatCount,ReverseComplement,OriginalSequence,RepeatInstances\n");
        
        if (filtered_repeats) {
            for (int i = 0; i < filtered_count; i++) {
                fprintf(detail_file, "%d,%d,%d,%s,%s,", 
                       filtered_repeats[i].position,
                       filtered_repeats[i].length,
                       filtered_repeats[i].count,
                       filtered_repeats[i].is_reverse ? "Yes" : "No",
                       filtered_repeats[i].orig_seq);
                
                // Add repeat examples
                if (filtered_repeats[i].num_examples > 0) {
                    for (int j = 0; j < filtered_repeats[i].num_examples; j++) {
                        fprintf(detail_file, "%s", filtered_repeats[i].repeat_examples[j]);
                        if (j < filtered_repeats[i].num_examples - 1) {
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
    
    // Free all allocated memory
    if (filtered_repeats) {
        free_repeat_patterns(filtered_repeats, filtered_count);
    }
    free(reference);
    free(query);
    
    return 0;
}
