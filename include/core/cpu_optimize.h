#ifndef CPU_OPTIMIZE_H
#define CPU_OPTIMIZE_H

// CPU 特性检测和优化设置
#ifdef __x86_64__
    #include <immintrin.h>
#endif

// AMD Ryzen 9 7940HX specific optimizations
#define CPU_RYZEN_7940HX
#define NUM_PHYSICAL_CORES 16
#define NUM_LOGICAL_CORES 32
#define L1_CACHE_SIZE 32768     // 32KB per core
#define L2_CACHE_SIZE 512000    // 512KB per core
#define L3_CACHE_SIZE 32768000  // 32MB shared

// Cache line size for Zen 4 architecture
#define CACHE_LINE_SIZE 64

// Prefetch settings
#define PREFETCH_DISTANCE (CACHE_LINE_SIZE * 2)
#define SOFTWARE_PREFETCH_HINT _MM_HINT_T0

// 内存对齐和预取宏
#define ALIGN_TO_CACHE __attribute__((aligned(CACHE_LINE_SIZE)))
#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

// 向量化和SIMD优化
#ifdef __AVX2__
    #define USE_AVX2
#endif

#ifdef __AVX512F__
    #define USE_AVX512
#endif

// 线程和并行计算设置
#ifdef _OPENMP
    #define DEFAULT_NUM_THREADS NUM_LOGICAL_CORES
    #define PARALLEL_THRESHOLD 1000  // 并行化的最小数据量阈值
#endif

// 内存预取函数
static inline void prefetch_read(const void* addr) {
    #ifdef __x86_64__
        _mm_prefetch((const char*)addr, SOFTWARE_PREFETCH_HINT);
    #endif
}

// 优化的内存分配函数
static inline void* aligned_alloc_cache(size_t size) {
    void* ptr = NULL;
    int ret = posix_memalign(&ptr, CACHE_LINE_SIZE, size);
    return (ret == 0) ? ptr : NULL;
}

// SIMD指令集检测
static inline int cpu_supports_avx2(void) {
    #ifdef __x86_64__
        unsigned int eax, ebx, ecx, edx;
        __get_cpuid(7, &eax, &ebx, &ecx, &edx);
        return (ebx & bit_AVX2) != 0;
    #else
        return 0;
    #endif
}

#endif // CPU_OPTIMIZE_H
