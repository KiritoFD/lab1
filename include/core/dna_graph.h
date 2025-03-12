#ifndef DNA_GRAPH_H
#define DNA_GRAPH_H

#include "dna_common.h"

// Structure to represent a node in our DNA graph
typedef struct GraphNode {
    int position;                // Position in sequence
    struct GraphEdge* edges;     // Outgoing edges
    int num_edges;               // Number of edges
    int capacity;                // Capacity of edges array
} GraphNode;

// Structure to represent an edge in our DNA graph
typedef struct GraphEdge {
    struct GraphNode* target;    // Target node
    int match_length;            // Length of matching subsequence
    double weight;               // Edge weight (based on match quality)
    int is_reverse;              // Whether this is a reverse complement match
} GraphEdge;

// Structure for our DNA graph
typedef struct {
    GraphNode* nodes;            // Array of nodes
    int num_nodes;               // Number of nodes
    int capacity;                // Capacity of nodes array
} DNAGraph;

// Function prototypes for graph operations
DNAGraph* build_dna_graph(const char* reference, int ref_len, const char* query, int query_len);
void add_edge(GraphNode* source, GraphNode* target, int match_length, int is_reverse);
RepeatPattern* find_repeats_in_graph(DNAGraph* graph, const char* reference, const char* query, int* num_repeats);
void free_dna_graph(DNAGraph* graph);

#endif // DNA_GRAPH_H
