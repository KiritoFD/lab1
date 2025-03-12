#include "../include/dna_graph.h"

// Build a directed acyclic graph representation of the DNA sequences
DNAGraph* build_dna_graph(const char* reference, int ref_len, const char* query, int query_len) {
    printf("Building DNA graph for pattern matching...\n");
    
    // Use all parameters to avoid warnings
    if (!reference || ref_len <= 0 || !query || query_len <= 0) {
        fprintf(stderr, "Invalid parameters for graph building\n");
        return NULL;
    }
    
    // Allocate and initialize the graph
    DNAGraph* graph = (DNAGraph*)malloc(sizeof(DNAGraph));
    if (!graph) {
        fprintf(stderr, "Failed to allocate memory for graph\n");
        return NULL;
    }
    
    // Create nodes for each position in the reference sequence
    int initial_capacity = ref_len;
    graph->nodes = (GraphNode*)malloc(initial_capacity * sizeof(GraphNode));
    if (!graph->nodes) {
        fprintf(stderr, "Failed to allocate memory for graph nodes\n");
        free(graph);
        return NULL;
    }
    
    graph->num_nodes = ref_len;
    graph->capacity = initial_capacity;
    
    // Initialize each node
    #pragma omp parallel for
    for (int i = 0; i < ref_len; i++) {
        graph->nodes[i].position = i;
        graph->nodes[i].edges = NULL;
        graph->nodes[i].num_edges = 0;
        graph->nodes[i].capacity = 0;
    }
    
    // Parameters for matching
    int min_length = 5 > (ref_len / 1000) ? 5 : (ref_len / 1000);
    
    // For large reference sequences, limit the positions we check
    int max_positions_to_check = 10000;
    int positions_step = ref_len > max_positions_to_check ? (ref_len / max_positions_to_check) : 1;
    
    printf("Building graph edges with min_length=%d...\n", min_length);
    
    // Build edges between nodes based on sequence matches
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < ref_len - min_length; i += positions_step) {
        // Extract segment from reference
        char* segment = (char*)malloc((min_length + 1) * sizeof(char));
        if (!segment) continue;
        
        strncpy(segment, reference + i, min_length);
        segment[min_length] = '\0';
        
        // Find matches in query
        int start_idx = 0;
        while (1) {
            char* found = strstr(query + start_idx, segment);
            if (!found) break;
            
            int match_pos = found - query;
            
            // Add an edge in the graph
            #pragma omp critical
            {
                add_edge(&graph->nodes[i], &graph->nodes[match_pos], min_length, 0);
            }
            
            start_idx = match_pos + 1;
            if (start_idx >= query_len - min_length) break;
        }
        
        // Check for reverse complement matches
        char* rev_comp = get_reverse_complement(segment, min_length);
        start_idx = 0;
        
        while (1) {
            char* found = strstr(query + start_idx, rev_comp);
            if (!found) break;
            
            int match_pos = found - query;
            
            // Add an edge for reverse complement match
            #pragma omp critical
            {
                add_edge(&graph->nodes[i], &graph->nodes[match_pos], min_length, 1);
            }
            
            start_idx = match_pos + 1;
            if (start_idx >= query_len - min_length) break;
        }
        
        free(segment);
        free(rev_comp);
    }
    
    printf("Graph construction complete: %d nodes created\n", graph->num_nodes);
    return graph;
}

// Add an edge between two nodes in the graph
void add_edge(GraphNode* source, GraphNode* target, int match_length, int is_reverse) {
    if (!source || !target) return; // Validate parameters
    
    // Allocate initial edges array if needed
    if (source->edges == NULL) {
        source->capacity = 5;  // Start with space for 5 edges
        source->edges = (GraphEdge*)malloc(source->capacity * sizeof(GraphEdge));
        if (!source->edges) return;
    }
    
    // Expand edges array if needed
    if (source->num_edges >= source->capacity) {
        source->capacity *= 2;
        GraphEdge* new_edges = (GraphEdge*)realloc(source->edges, source->capacity * sizeof(GraphEdge));
        if (!new_edges) return;
        source->edges = new_edges;
    }
    
    // Add the new edge with parameters
    source->edges[source->num_edges].target = target;
    source->edges[source->num_edges].match_length = match_length;
    source->edges[source->num_edges].weight = match_length + (is_reverse ? 0.5 : 0);  // Weight based on match length
    source->edges[source->num_edges].is_reverse = is_reverse;
    source->num_edges++;
}

// Find repeats by traversing paths in the DNA graph
RepeatPattern* find_repeats_in_graph(DNAGraph* graph, const char* reference, const char* query, int* num_repeats) {
    if (!graph || !reference || !query || !num_repeats) {
        if (num_repeats) *num_repeats = 0;
        return NULL;
    }
    
    printf("Finding repeats using graph traversal...\n");
    
    // Just return a placeholder for now until fully implemented
    *num_repeats = 0;
    return NULL;
}

// Free memory used by the DNA graph
void free_dna_graph(DNAGraph* graph) {
    if (!graph) return;
    
    if (graph->nodes) {
        for (int i = 0; i < graph->num_nodes; i++) {
            if (graph->nodes[i].edges) {
                free(graph->nodes[i].edges);
            }
        }
        free(graph->nodes);
    }
    
    free(graph);
}
