#include "../include/core/dna_common.h"
#include "../include/core/dna_io.h"
#include "../include/core/dna_traditional.h"
#include "../include/core/dna_graph.h"
#include <sys/stat.h>
#include <time.h>

// Default input file paths
#define DEFAULT_REFERENCE_FILE "reference.txt"
#define DEFAULT_QUERY_FILE "query.txt"
#define OUTPUT_DIR "output"

// Function to check if directory exists, create if not
void ensure_output_directory() {
    struct stat st = {0};
    if (stat(OUTPUT_DIR, &st) == -1) {
        #ifdef _WIN32
            mkdir(OUTPUT_DIR);
        #else
            mkdir(OUTPUT_DIR, 0700);
        #endif
        printf("Created output directory: %s\n", OUTPUT_DIR);
    }
}

// Function to get next available file number
int get_next_file_number() {
    int file_num = 1;
    char filepath[100];
    FILE *file;
    
    // Find the next available file number
    while (1) {
        snprintf(filepath, sizeof(filepath), "%s/result_%d.txt", OUTPUT_DIR, file_num);
        file = fopen(filepath, "r");
        if (file == NULL) {
            break;  // File doesn't exist, use this number
        }
        fclose(file);
        file_num++;
    }
    
    return file_num;
}

// The main function
int main(int argc, char *argv[]) {
    // Ensure output directory exists
    ensure_output_directory();
    
    // Get next available file number
    int file_num = get_next_file_number();
    char output_filepath[100];
    snprintf(output_filepath, sizeof(output_filepath), "%s/result_%d.txt", OUTPUT_DIR, file_num);
    
    // Create output file
    FILE *output_file = fopen(output_filepath, "w");
    if (!output_file) {
        fprintf(stderr, "Failed to create output file: %s\n", output_filepath);
        return EXIT_FAILURE;
    }
    
    char* reference_file;
    char* query_file;
    
    // Check if file paths are provided as arguments, otherwise use defaults
    if (argc == 3) {
        reference_file = argv[1];
        query_file = argv[2];
        printf("Using provided file paths:\n");
        fprintf(output_file, "Using provided file paths:\n");
    } else {
        reference_file = DEFAULT_REFERENCE_FILE;
        query_file = DEFAULT_QUERY_FILE;
        printf("Using default file paths:\n");
        fprintf(output_file, "Using default file paths:\n");
    }
    
    printf("Reference file: %s\n", reference_file);
    printf("Query file: %s\n", query_file);
    fprintf(output_file, "Reference file: %s\n", reference_file);
    fprintf(output_file, "Query file: %s\n", query_file);
    
    printf("DNA Repeat Finder\n");
    printf("Version: 2.0 with DAG-based approach\n\n");
    fprintf(output_file, "DNA Repeat Finder\n");
    fprintf(output_file, "Version: 2.0 with DAG-based approach\n\n");
    
    // Thread handling (if OpenMP is enabled)
    #ifdef _OPENMP
    int max_threads = omp_get_max_threads();
    int thread_count = max_threads; // Default to maximum
    
    // Automatically use maximum threads without asking
    printf("Using %d threads for parallel processing\n", thread_count);
    fprintf(output_file, "Using %d threads for parallel processing\n", thread_count);
    omp_set_num_threads(thread_count);
    #endif
    
    char* reference = NULL;
    char* query = NULL;
    int ref_len = 0;
    int query_len = 0;
    
    // Read sequences from the provided files
    printf("Reading reference sequence from %s...\n", reference_file);
    fprintf(output_file, "Reading reference sequence from %s...\n", reference_file);
    reference = read_sequence_from_file(reference_file, &ref_len);
    if (!reference) {
        fprintf(stderr, "Failed to read reference file: %s\n", reference_file);
        fprintf(output_file, "Failed to read reference file: %s\n", reference_file);
        fclose(output_file);
        return EXIT_FAILURE;
    }
    
    printf("Reading query sequence from %s...\n", query_file);
    fprintf(output_file, "Reading query sequence from %s...\n", query_file);
    query = read_sequence_from_file(query_file, &query_len);
    if (!query) {
        fprintf(stderr, "Failed to read query file: %s\n", query_file);
        fprintf(output_file, "Failed to read query file: %s\n", query_file);
        free(reference);
        fclose(output_file);
        return EXIT_FAILURE;
    }
    
    printf("Successfully loaded sequences. Reference length: %d, Query length: %d\n", ref_len, query_len);
    fprintf(output_file, "Successfully loaded sequences. Reference length: %d, Query length: %d\n", ref_len, query_len);
    
    // Record start time for graph approach
    clock_t start_time = clock();
    
    // Build DNA graph and find repeats using graph-based approach
    printf("\n--- Using DAG-based approach ---\n");
    fprintf(output_file, "\n--- Using DAG-based approach ---\n");
    DNAGraph* dna_graph = build_dna_graph(reference, ref_len, query, query_len);
    
    int num_graph_repeats = 0;
    RepeatPattern* graph_repeats = find_repeats_in_graph(dna_graph, reference, query, &num_graph_repeats);
    
    // Free graph memory
    free_dna_graph(dna_graph);
    
    clock_t graph_end_time = clock();
    double graph_time = ((double)(graph_end_time - start_time) * 1000.0) / CLOCKS_PER_SEC;
    
    // Filter nested repeats for graph-based approach
    int filtered_graph_count = 0;
    RepeatPattern* filtered_graph_repeats = NULL;
    
    if (graph_repeats) {
        // First get sequence information
        int seq_count = 0;
        RepeatPattern* repeats_with_seq = get_repeat_sequences(graph_repeats, num_graph_repeats, 
                                                              reference, query, ref_len, query_len, &seq_count);
        
        // Then filter nested repeats
        filtered_graph_repeats = filter_nested_repeats(repeats_with_seq, seq_count, 1, &filtered_graph_count);
    }
    
    // Display graph-based results
    printf("\nGraph-based approach found %d unique repeat patterns\n", filtered_graph_count);
    printf("Graph processing time: %.2f milliseconds\n", graph_time);
    fprintf(output_file, "\nGraph-based approach found %d unique repeat patterns\n", filtered_graph_count);
    fprintf(output_file, "Graph processing time: %.2f milliseconds\n", graph_time);
    
    // Save graph-based results to console and file
    if (filtered_graph_repeats) {
        for (int i = 0; i < filtered_graph_count; i++) {
            printf("Repeat Pattern %d: Position: %d, Length: %d, Count: %d, Is Reverse: %d\n", 
                  i+1, filtered_graph_repeats[i].position, filtered_graph_repeats[i].length, 
                  filtered_graph_repeats[i].count, filtered_graph_repeats[i].is_reverse);
            fprintf(output_file, "Repeat Pattern %d: Position: %d, Length: %d, Count: %d, Is Reverse: %d\n", 
                   i+1, filtered_graph_repeats[i].position, filtered_graph_repeats[i].length, 
                   filtered_graph_repeats[i].count, filtered_graph_repeats[i].is_reverse);
            
            // Print original sequence if available
            if (filtered_graph_repeats[i].orig_seq) {
                printf("  Sequence: %s\n", filtered_graph_repeats[i].orig_seq);
                fprintf(output_file, "  Sequence: %s\n", filtered_graph_repeats[i].orig_seq);
            }
            
            // Print repeat examples if available
            for (int j = 0; j < filtered_graph_repeats[i].num_examples && j < 3; j++) {
                if (filtered_graph_repeats[i].repeat_examples[j]) {
                    printf("  Example %d: %s\n", j+1, filtered_graph_repeats[i].repeat_examples[j]);
                    fprintf(output_file, "  Example %d: %s\n", j+1, filtered_graph_repeats[i].repeat_examples[j]);
                }
            }
        }
        free_repeat_patterns(filtered_graph_repeats, filtered_graph_count);
    }
    
    printf("\nResults have been saved to: %s\n", output_filepath);
    
    // Close output file
    fclose(output_file);
    
    // Free memory
    free(reference);
    free(query);
    
    return 0;
}
