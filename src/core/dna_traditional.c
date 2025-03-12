#include "../include/core/dna_traditional.h"
#include "../include/core/cpu_optimize.h"

// Build similarity matrix between reference and query - optimized with parallel processing and AVX2
int** build_similarity_matrix(const char* reference, int ref_len, const char* query, int query_len) {
    int** matrix = (int**)aligned_alloc(CACHE_LINE_SIZE, ref_len * sizeof(int*));
    if (UNLIKELY(!matrix)) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    // Determine optimal thread count for this workload
    int thread_count = get_optimal_thread_count((size_t)ref_len * query_len * sizeof(int));
    omp_set_num_threads(thread_count);
    
    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < ref_len; i++) {
        // Cache line aligned allocation for better memory access
        matrix[i] = (int*)aligned_alloc(CACHE_LINE_SIZE, ((query_len * sizeof(int) + CACHE_LINE_SIZE - 1) 
                                                         & ~(CACHE_LINE_SIZE - 1)));
        if (UNLIKELY(!matrix[i])) {
            fprintf(stderr, "Memory allocation failed in thread %d\n", omp_get_thread_num());
            // Cannot safely free memory here due to parallel context
            exit(EXIT_FAILURE);
        }
        
        // Prefetch next rows to improve cache utilization
        if (i + 2 < ref_len) {
            PREFETCH_READ(&reference[i + 2]);
        }

        #ifdef __AVX2__
        // Process 16 elements at a time using AVX2
        for (int j = 0; j < query_len - 15; j += 16) {
            // Load 16 characters from query
            __m256i q_chars1 = _mm256_loadu_si256((__m256i*)&query[j]);
            __m256i q_chars2 = _mm256_loadu_si256((__m256i*)&query[j+8]);
            
            // Compare with current reference character
            __m256i ref_char1 = _mm256_set1_epi8(reference[i]);
            __m256i ref_char2 = _mm256_set1_epi8(reference[i]);
            
            __m256i match1 = _mm256_cmpeq_epi8(q_chars1, ref_char1);
            __m256i match2 = _mm256_cmpeq_epi8(q_chars2, ref_char2);
            
            // Remove unused variables
            // Convert to int (-1 for mismatch, 1 for match)
            
            // Process first 8 matches
            int mask1 = _mm256_movemask_epi8(match1);
            for (int k = 0; k < 8; k++) {
                matrix[i][j+k] = (mask1 & (1 << k)) ? 1 : -1;
            }
            
            // Process next 8 matches
            int mask2 = _mm256_movemask_epi8(match2);
            for (int k = 0; k < 8; k++) {
                matrix[i][j+8+k] = (mask2 & (1 << k)) ? 1 : -1;
            }
        }
        
        // Handle remaining elements
        for (int j = (query_len & ~15); j < query_len; j++) {
            matrix[i][j] = (reference[i] == query[j]) ? 1 : -1;
        }
        #else
        // Fallback for non-AVX2 systems
        for (int j = 0; j < query_len; j++) {
            matrix[i][j] = (reference[i] == query[j]) ? 1 : -1;
        }
        #endif
    }
    
    return matrix;
}

// Find repeats in the query sequence compared to reference - optimized for R9 7940HX
RepeatPattern* find_repeats(const char* reference, int ref_len, const char* query, int query_len, int* num_repeats) {
    // Dynamic parameters based on sequence length with optimized defaults for performance
    int min_length = 5 > (ref_len / 1000) ? 5 : (ref_len / 1000);
    int max_length = (ref_len / 10) < 120 ? (ref_len / 10) : 120;
    int step = 1 > (min_length / 5) ? 1 : (min_length / 5);
    
    printf("Using parameters: min_length=%d, max_length=%d, step=%d\n", 
           min_length, max_length, step);
    
    // Limit search space for very large sequences
    int max_positions_to_check = 10000;
    int positions_step = ref_len > max_positions_to_check ? (ref_len / max_positions_to_check) : 1;
    
    // Allocate initial memory for repeats (aligned for better cache performance)
    int capacity = 1000;
    RepeatPattern* repeats = (RepeatPattern*)aligned_alloc(CACHE_LINE_SIZE, 
                                            capacity * sizeof(RepeatPattern));
    if (UNLIKELY(!repeats)) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    int repeat_count = 0;
    
    // Progress reporting
    int progress = 0;
    int total_progress = (ref_len / positions_step);
    int last_reported_progress = -1;
    
    // Optimize thread count based on workload
    int thread_count = get_optimal_thread_count((size_t)ref_len * query_len);
    omp_set_num_threads(thread_count);
    printf("Using %d threads for repeat finding\n", thread_count);

    // For thread-safe updates to repeats array
    omp_lock_t repeat_lock;
    omp_init_lock(&repeat_lock);
    
    #pragma omp parallel
    {
        // Thread-local repeats for better performance
        int local_capacity = 100;
        RepeatPattern* local_repeats = (RepeatPattern*)malloc(local_capacity * sizeof(RepeatPattern));
        int local_count = 0;
        
        #pragma omp for schedule(dynamic, 16) 
        for (int pos = 0; pos < ref_len - min_length; pos += positions_step) {
            // Report progress (only from thread 0)
            if (omp_get_thread_num() == 0) {
                progress++;
                int current_progress = (progress * 100) / total_progress;
                if (current_progress % 10 == 0 && current_progress != last_reported_progress) {
                    printf("Finding repeats: %d%% complete\n", current_progress);
                    last_reported_progress = current_progress;
                }
            }
            
            // Check different segment lengths
            for (int length = min_length; length < (max_length < (ref_len - pos) ? max_length : (ref_len - pos)); length += step) {
                // Extract segment from reference with cache alignment
                char* segment = allocate_dna_sequence(length + 1);
                strncpy(segment, reference + pos, length);
                segment[length] = '\0';
                
                // Pre-compute hash for faster string matching
                unsigned long segment_hash = 0;
                for (int i = 0; i < length; i++) {
                    segment_hash = segment_hash * 31 + segment[i];
                }
                
                // Check for forward repeats with optimized scanning
                int start_idx = 0;
                while (1) {
                    // Find next occurrence of segment
                    char* found = strstr(query + start_idx, segment);
                    if (!found) break;
                    
                    int next_idx = found - query;
                    int current_pos = next_idx + length;
                    int consecutive_count = 0;
                    
                    // Prefetch ahead for better cache performance
                    PREFETCH_READ(query + current_pos);
                    
                    // Check for consecutive repeats using vectorized comparison when available
                    while (current_pos + length <= query_len) {
                        int is_repeat;
                        
                        #ifdef __AVX2__
                        if (length >= 32) {
                            is_repeat = vectorized_dna_compare(query + current_pos, segment, length);
                        } else {
                        #endif
                            is_repeat = 1;
                            for (int k = 0; k < length; k++) {
                                if (query[current_pos + k] != segment[k]) {
                                    is_repeat = 0;
                                    break;
                                }
                            }
                        #ifdef __AVX2__
                        }
                        #endif
                        
                        if (is_repeat) {
                            consecutive_count++;
                            current_pos += length;
                            // Prefetch next position for better performance
                            if (current_pos + length <= query_len)
                                PREFETCH_READ(query + current_pos);
                        } else {
                            break;
                        }
                    }
                    
                    if (consecutive_count > 0) {
                        // Add to local repeats array
                        if (local_count >= local_capacity) {
                            local_capacity *= 2;
                            RepeatPattern* new_local = (RepeatPattern*)realloc(local_repeats, 
                                                     local_capacity * sizeof(RepeatPattern));
                            if (UNLIKELY(!new_local)) {
                                fprintf(stderr, "Memory reallocation failed\n");
                                free(segment);
                                free(local_repeats);
                                exit(EXIT_FAILURE);
                            }
                            local_repeats = new_local;
                        }
                        
                        local_repeats[local_count].position = pos;
                        local_repeats[local_count].length = length;
                        local_repeats[local_count].count = consecutive_count;
                        local_repeats[local_count].is_reverse = 0;
                        local_repeats[local_count].orig_seq = strdup(segment);
                        local_repeats[local_count].repeat_examples = NULL;
                        local_repeats[local_count].num_examples = 0;
                        local_count++;
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
                        int is_repeat;
                        
                        #ifdef __AVX2__
                        if (length >= 32) {
                            is_repeat = vectorized_dna_compare(query + current_pos, rev_comp, length);
                        } else {
                        #endif
                            is_repeat = 1;
                            for (int k = 0; k < length; k++) {
                                if (query[current_pos + k] != rev_comp[k]) {
                                    is_repeat = 0;
                                    break;
                                }
                            }
                        #ifdef __AVX2__
                        }
                        #endif
                        
                        if (is_repeat) {
                            consecutive_count++;
                            current_pos += length;
                        } else {
                            break;
                        }
                    }
                    
                    // Add to local repeats array
                    if (local_count >= local_capacity) {
                        local_capacity *= 2;
                        RepeatPattern* new_local = (RepeatPattern*)realloc(local_repeats, 
                                                 local_capacity * sizeof(RepeatPattern));
                        if (UNLIKELY(!new_local)) {
                            fprintf(stderr, "Memory reallocation failed\n");
                            free(segment);
                            free(rev_comp);
                            free(local_repeats);
                            exit(EXIT_FAILURE);
                        }
                        local_repeats = new_local;
                    }
                    
                    local_repeats[local_count].position = pos;
                    local_repeats[local_count].length = length;
                    local_repeats[local_count].count = consecutive_count > 0 ? consecutive_count : 1;
                    local_repeats[local_count].is_reverse = 1;
                    local_repeats[local_count].orig_seq = strdup(segment);
                    local_repeats[local_count].repeat_examples = NULL;
                    local_repeats[local_count].num_examples = 0;
                    local_count++;
                    
                    start_idx = next_idx + 1;
                    if (start_idx >= query_len) break;
                }
                
                free(segment);
                free(rev_comp);
            }
        }
        
        // Merge thread-local results into global array
        omp_set_lock(&repeat_lock);
        if (repeat_count + local_count > capacity) {
            // Need to resize global array
            while(repeat_count + local_count > capacity) capacity *= 2;
            RepeatPattern* new_repeats = (RepeatPattern*)realloc(repeats, capacity * sizeof(RepeatPattern));
            if (UNLIKELY(!new_repeats)) {
                fprintf(stderr, "Memory reallocation failed during merge\n");
                free(local_repeats);
                omp_unset_lock(&repeat_lock);
                exit(EXIT_FAILURE);
            }
            repeats = new_repeats;
        }
        
        // Copy local results to global array
        memcpy(repeats + repeat_count, local_repeats, local_count * sizeof(RepeatPattern));
        repeat_count += local_count;
        omp_unset_lock(&repeat_lock);
        
        // Free thread-local array
        free(local_repeats);
    }
    
    omp_destroy_lock(&repeat_lock);
    
    printf("Found %d repeat patterns\n", repeat_count);
    *num_repeats = repeat_count;
    
    // If no repeats found, free memory and return NULL
    if (repeat_count == 0) {
        free(repeats);
        return NULL;
    }
    
    return repeats;
}

// Add sequence information to repeat patterns - vectorized and parallelized version
RepeatPattern* get_repeat_sequences(RepeatPattern* repeats, int num_repeats, const char* reference, const char* query, int ref_len, int query_len, int* new_count) {
    // Remove unused parameter warning with attribute
    (void)ref_len; // Explicitly mark as unused
    
    *new_count = num_repeats;
    
    // Set optimal thread count based on workload
    int thread_count = get_optimal_thread_count(num_repeats * 100);
    omp_set_num_threads(thread_count);
    
    #pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < num_repeats; i++) {
        int pos = repeats[i].position;
        int length = repeats[i].length;
        int is_reverse = repeats[i].is_reverse;
        
        // Get the original segment from reference
        char* segment = allocate_dna_sequence(length + 1);
        strncpy(segment, reference + pos, length);
        segment[length] = '\0';
        
        // For reverse complement
        char* search_seq;
        if (is_reverse) {
            search_seq = get_reverse_complement(segment, length);
        } else {
            search_seq = strdup(segment);
        }
        
        // Find repeat instances in query with optimized search
        int start_idx = 0;
        int instances_count = 0;
        int instances_capacity = 16; // Start with reasonable capacity
        char** instances = (char**)malloc(instances_capacity * sizeof(char*));
        
        while (1) {
            // Optimize string search
            char* found = strstr(query + start_idx, search_seq);
            if (!found) break;
            
            int next_idx = found - query;
            
            // Resize instances array if needed
            if (instances_count >= instances_capacity) {
                instances_capacity *= 2;
                char** new_instances = (char**)realloc(instances, instances_capacity * sizeof(char*));
                if (UNLIKELY(!new_instances)) {
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
            
            // Store the found instance
            instances[instances_count] = allocate_dna_sequence(length + 1);
            strncpy(instances[instances_count], query + next_idx, length);
            instances[instances_count][length] = '\0';
            instances_count++;
            
            // Advance search position
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
        RepeatPattern* temp = (RepeatPattern*)aligned_alloc(CACHE_LINE_SIZE, 
                                           num_repeats * sizeof(RepeatPattern));
        if (!temp) {
            fprintf(stderr, "Memory allocation failed in filtering\n");
            *filtered_count = num_repeats;
            return repeats;
        }
        
        // First pass - count valid entries
        #pragma omp parallel for reduction(+:count)
        for (int i = 0; i < num_repeats; i++) {
            if (repeats[i].num_examples > 0) {
                count++;
            }
        }
        
        // Second pass - copy valid entries
        int idx = 0;
        for (int i = 0; i < num_repeats; i++) {
            if (repeats[i].num_examples > 0) {
                temp[idx++] = repeats[i];
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
    
    PositionKey* positions = (PositionKey*)aligned_alloc(CACHE_LINE_SIZE, 
                                          filtered_size * sizeof(PositionKey));
    if (!positions) {
        fprintf(stderr, "Memory allocation failed in position tracking\n");
        if (filtered != repeats) free(filtered);
        *filtered_count = num_repeats;
        return repeats;
    }
    
    // Initialize position keys
    #pragma omp parallel for
    for (int i = 0; i < filtered_size; i++) {
        positions[i].position = filtered[i].position;
        positions[i].is_reverse = filtered[i].is_reverse;
        positions[i].used = 0;
    }
    
    // Allocate result array (aligned for better performance)
    RepeatPattern* result = (RepeatPattern*)aligned_alloc(CACHE_LINE_SIZE, 
                                          filtered_size * sizeof(RepeatPattern));
    if (!result) {
        fprintf(stderr, "Memory allocation failed for results\n");
        free(positions);
        if (filtered != repeats) free(filtered);
        *filtered_count = num_repeats;
        return repeats;
    }
    
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

// Free memory used by similarity matrix - optimized for large matrices
void free_matrix(int** matrix, int rows) {
    if (!matrix) return;
    
    // Use parallel processing for large matrices
    if (rows > 1000) {
        #pragma omp parallel for
        for (int i = 0; i < rows; i++) {
            if (matrix[i]) free(matrix[i]);
        }
    } else {
        // Sequential for smaller matrices
        for (int i = 0; i < rows; i++) {
            if (matrix[i]) free(matrix[i]);
        }
    }
    
    free(matrix);
}
