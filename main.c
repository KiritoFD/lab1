#include "../include/dna_common.h"
#include "../include/dna_io.h"
#include "../include/dna_traditional.h"
#include "../include/dna_graph.h"

// Default input file paths
#define DEFAULT_REFERENCE_FILE "reference.txt"
#define DEFAULT_QUERY_FILE "query.txt"

// Main function
int main(int argc, char *argv[]) {
    char* reference_file;
    char* query_file;
    
    // Check if file paths are provided as arguments, otherwise use defaults
    if (argc == 3) {
        reference_file = argv[1];
        query_file = argv[2];
        printf("Using provided file paths:\n");
    } else {
        reference_file = DEFAULT_REFERENCE_FILE;
        query_file = DEFAULT_QUERY_FILE;
        printf("Using default file paths:\n");
    }
    
    printf("Reference file: %s\n", reference_file);
    printf("Query file: %s\n", query_file);
    
    printf("DNA Repeat Finder\n");
    printf("Version: 2.0 with DAG-based approach\n\n");
    
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
    
    // Record start time for graph approach
    clock_t start_time = clock();
    
    // Build DNA graph and find repeats using graph-based approach
    printf("\n--- Using DAG-based approach ---\n");
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
    
    // Save graph-based results
    if (filtered_graph_repeats) {
        save_results(filtered_graph_repeats, filtered_graph_count, num_graph_repeats);
        free_repeat_patterns(filtered_graph_repeats, filtered_graph_count);
    }
    
    // Free memory
    free(reference);
    free(query);
    
    return 0;
}
