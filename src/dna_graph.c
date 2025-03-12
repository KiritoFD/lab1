#include "../include/dna_graph.h"

// Build a directed acyclic graph representation of the DNA sequences
DNAGraph* build_dna_graph(const char* reference, int ref_len, const char* query, int query_len) {
    printf("Building DNA graph for pattern matching...\n");
    
    // Allocate and initialize the graph
    DNAGraph* graph = (DNAGraph*)malloc(sizeof(DNAGraph));
    if (!graph) {
        fprintf(stderr, "Failed to allocate memory for graph\n");
        exit(EXIT_FAILURE);
    }
    
    // Create nodes for each position in the reference sequence
    int initial_capacity = ref_len;
    graph->nodes = (GraphNode*)malloc(initial_capacity * sizeof(GraphNode));
    if (!graph->nodes) {
        fprintf(stderr, "Failed to allocate memory for graph nodes\n");
        free(graph);
        exit(EXIT_FAILURE);
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
    
    // Add the new edge
    source->edges[source->num_edges].target = target;
    source->edges[source->num_edges].match_length = match_length;
    source->edges[source->num_edges].weight = match_length;  // Weight based on match length
    source->edges[source->num_edges].is_reverse = is_reverse;
    source->num_edges++;
}

// Find repeats by traversing paths in the DNA graph
RepeatPattern* find_repeats_in_graph(DNAGraph* graph, const char* reference, const char* query, int* num_repeats) {
    printf("Finding repeats using graph traversal...\n");
    
    // Allocate initial memory for repeats
    int capacity = 1000;
    RepeatPattern* repeats = (RepeatPattern*)aligned_alloc(CACHE_LINE_SIZE,
                                            capacity * sizeof(RepeatPattern));
    if (!repeats) {
        fprintf(stderr, "Memory allocation failed for repeats\n");
        *num_repeats = 0;
        return NULL;
    }
    
    int repeat_count = 0;
    
    // Track visited nodes to avoid cycles (though our graph should be acyclic)
    int* visited = (int*)calloc(graph->num_nodes, sizeof(int));
    if (!visited) {
        free(repeats);
        *num_repeats = 0;
        return NULL;
    }
    
    // For each node, find paths that represent repeats
    for (int i = 0; i < graph->num_nodes; i++) {
        GraphNode* node = &graph->nodes[i];
        
        // Skip nodes with no outgoing edges
        if (node->num_edges == 0) continue;
        
        // For each edge, check if it forms a repeat pattern
        for (int j = 0; j < node->num_edges; j++) {
            GraphEdge* edge = &node->edges[j];
            int target_pos = edge->target->position;
            int match_length = edge->match_length;
            int is_reverse = edge->is_reverse;
            
            // Extract the original sequence
            char* orig_seq = (char*)malloc((match_length + 1) * sizeof(char));
            strncpy(orig_seq, reference + node->position, match_length);
            orig_seq[match_length] = '\0';
            
            // Count consecutive repeats
            int consecutive_count = 0;
            int current_pos = target_pos + match_length;
            
            // For forward matches
            if (!is_reverse) {
                size_t query_len = strlen(query);
                while ((size_t)(current_pos + match_length) <= query_len) {
                    int is_match = 1;
                    for (int k = 0; k < match_length; k++) {
                        if (query[current_pos + k] != orig_seq[k]) {
                            is_match = 0;
                            break;
                        }
                    }
                    
                    if (is_match) {
                        consecutive_count++;
                        current_pos += match_length;
                    } else {
                        break;
                    }
                }
            }
            // For reverse complement matches
            else {
                char* rev_comp = get_reverse_complement(orig_seq, match_length);
                size_t query_len = strlen(query);
                
                while ((size_t)(current_pos + match_length) <= query_len) {
                    int is_match = 1;
                    for (int k = 0; k < match_length; k++) {
                        if (query[current_pos + k] != rev_comp[k]) {
                            is_match = 0;
                            break;
                        }
                    }
                    
                    if (is_match) {
                        consecutive_count++;
                        current_pos += match_length;
                    } else {
                        break;
                    }
                }
                
                free(rev_comp);
            }
            
            // Add to repeats if we found consecutive matches or any reverse complement match
            if (consecutive_count > 0 || is_reverse) {
                // Resize repeats array if needed
                if (repeat_count >= capacity) {
                    capacity *= 2;
                    RepeatPattern* new_repeats = (RepeatPattern*)realloc(repeats, capacity * sizeof(RepeatPattern));
                    if (!new_repeats) {
                        free(orig_seq);
                        free(visited);
                        free_repeat_patterns(repeats, repeat_count);
                        *num_repeats = 0;
                        return NULL;
                    }
                    repeats = new_repeats;
                }
                
                // Add the repeat pattern
                repeats[repeat_count].position = node->position;
                repeats[repeat_count].length = match_length;
                repeats[repeat_count].count = consecutive_count > 0 ? consecutive_count : 1;
                repeats[repeat_count].is_reverse = is_reverse;
                repeats[repeat_count].orig_seq = orig_seq;
                repeats[repeat_count].repeat_examples = NULL;
                repeats[repeat_count].num_examples = 0;
                
                repeat_count++;
            } else {
                free(orig_seq);
            }
        }
    }
    
    free(visited);
    
    printf("Found %d repeat patterns using graph method\n", repeat_count);
    *num_repeats = repeat_count;
    
    if (repeat_count == 0) {
        free(repeats);
        return NULL;
    }
    
    return repeats;
}

// Free memory used by the DNA graph
void free_dna_graph(DNAGraph* graph) {
    if (!graph) return;
    
    for (int i = 0; i < graph->num_nodes; i++) {
        if (graph->nodes[i].edges) {
            free(graph->nodes[i].edges);
        }
    }
    
    free(graph->nodes);
    free(graph);
}
