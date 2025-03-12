#ifndef CPU_OPTIMIZE_H
#define CPU_OPTIMIZE_H

#include <immintrin.h> // For AVX2, FMA instructions

// CPU-specific optimization settings for AMD Ryzen 9 7940HX
#define R9_7940HX_L1_CACHE_SIZE (32 * 1024)      // 32KB L1 cache 
#define R9_7940HX_L2_CACHE_SIZE (1024 * 1024)    // 1MB L2 cache
#define R9_7940HX_L3_CACHE_SIZE (16 * 1024 * 1024) // 16MB L3 cache
#define MAX_THREADS 16                           // R9 7940HX has 16 threads

// Function to set optimal thread count based on workload size
static inline int get_optimal_thread_count(size_t data_size) {
    if (data_size < R9_7940HX_L1_CACHE_SIZE) {
        return 1; // Too small for parallelization
    } else if (data_size < R9_7940HX_L2_CACHE_SIZE * 4) {
        return 4; // Small workload
    } else if (data_size < R9_7940HX_L3_CACHE_SIZE) {
        return 8; // Medium workload
    } else {
        return MAX_THREADS; // Large workload, use all cores
    }
}

// Cache-friendly memory alignment
#define CACHE_LINE_SIZE 64
#define ALIGN_TO_CACHE __attribute__((aligned(CACHE_LINE_SIZE)))

// Prefetch hints for optimal memory access patterns
#define PREFETCH_READ(addr) _mm_prefetch((const char*)(addr), _MM_HINT_T0)
#define PREFETCH_WRITE(addr) _mm_prefetch((const char*)(addr), _MM_HINT_ET0)

// Tell the compiler which branch is more likely to be taken
#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

// Force function inlining for performance-critical functions
#define FORCE_INLINE __attribute__((always_inline)) inline

// Vectorization helpers for Zen 4 architecture
#ifdef __AVX2__
// Function for vectorized nucleotide comparison
// This could be used to accelerate DNA sequence matching
FORCE_INLINE int vectorized_dna_compare(const char* seq1, const char* seq2, int length) {
    int i = 0;
    
    // Use AVX2 for 32-byte chunks
    for (; i + 32 <= length; i += 32) {
        __m256i a = _mm256_loadu_si256((const __m256i*)(seq1 + i));
        __m256i b = _mm256_loadu_si256((const __m256i*)(seq2 + i));
        
        __m256i cmp = _mm256_cmpeq_epi8(a, b);
        unsigned int mask = (unsigned int)_mm256_movemask_epi8(cmp);
        
        if (mask != 0xFFFFFFFFU)
            return 0;
    }
    
    // Handle remaining bytes
    for (; i < length; i++) {
        if (seq1[i] != seq2[i])
            return 0;
    }
    
    return 1;
}
#endif // __AVX2__

#endif // CPU_OPTIMIZE_H
