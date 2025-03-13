/*
针对 AMD Ryzen 9 7940HX 的优化编译命令:
gcc -march=znver4 -mtune=znver4 -O3 -ffast-math -flto -fuse-linker-plugin \
    -fprefetch-loop-arrays -funroll-loops -fomit-frame-pointer -mavx2 -mfma \
    -pthread -fopenmp \
    -DCPU_RYZEN_7940HX -DNUM_CORES=16 -DNUM_THREADS=32 \
    -DL1_CACHE_SIZE=32768 -DL2_CACHE_SIZE=512000 -DL3_CACHE_SIZE=32768000 \
    dna_repeat_finder_new.c -o dna_repeat_finder_new
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include <immintrin.h>

#define NUM_PHYSICAL_CORES 16
#define NUM_LOGICAL_CORES 32
#define L1_CACHE_SIZE 32768
#define L2_CACHE_SIZE 512000
#define L3_CACHE_SIZE 32768000
#define CACHE_LINE_SIZE 64
#define MAX_LENGTH 120
#define MIN_LENGTH 10
#define MAX_REPEATS 1000
#define HASH_TABLE_SIZE 32768  // 增大哈希表以减少冲突
#define BLOCK_SIZE 64  // 缓存行大小对齐的块大小

// 数据结构缓存对齐
#define ALIGN_TO_CACHE __attribute__((aligned(CACHE_LINE_SIZE)))

// 重复片段的数据结构
typedef struct ALIGN_TO_CACHE {
    int position;
    int length;
    int repeat_count;
    bool is_reverse;
    char* original_sequence;
    int query_position;
} RepeatPattern;

// 优化的哈希表节点结构
typedef struct HashNode {
    char* key;
    int* positions;
    int count;
    int capacity;
    struct HashNode* next;
} HashNode ALIGN_TO_CACHE;

// 优化的哈希表结构
typedef struct {
    HashNode** buckets;
    int size;
} HashMap ALIGN_TO_CACHE;

// 线程本地存储
static __thread char* thread_local_buffer = NULL;
static __thread size_t thread_local_buffer_size = 0;

// SIMD优化的字符串比较
static inline int avx2_strcmp(const char* str1, const char* str2, size_t len) {
    size_t i = 0;
    
    // 使用AVX2处理32字节块
    for (; i + 32 <= len; i += 32) {
        __m256i v1 = _mm256_loadu_si256((__m256i*)(str1 + i));
        __m256i v2 = _mm256_loadu_si256((__m256i*)(str2 + i));
        __m256i vcmp = _mm256_cmpeq_epi8(v1, v2);
        int mask = _mm256_movemask_epi8(vcmp);
        if (mask != 0xFFFFFFFF) {
            return 0;
        }
    }
    
    // 处理剩余字节
    return memcmp(str1 + i, str2 + i, len - i) == 0;
}

// 优化的哈希函数
static inline unsigned int optimized_hash(const char* str, int size) {
    unsigned int hash = 5381;
    size_t len = strlen(str);
    size_t i = 0;
    
    // 使用SIMD优化哈希计算
    for (; i + 32 <= len; i += 32) {
        __m256i block = _mm256_loadu_si256((__m256i*)(str + i));
        __m256i multiplier = _mm256_set1_epi32(33);
        __m256i hash_vec = _mm256_set1_epi32(hash);
        
        // 计算哈希值
        hash_vec = _mm256_mullo_epi32(hash_vec, multiplier);
        hash_vec = _mm256_add_epi32(hash_vec, 
            _mm256_cvtepu8_epi32(_mm256_extracti128_si256(block, 0)));
        
        // 合并结果
        int* hash_array = (int*)&hash_vec;
        for (int j = 0; j < 8; j++) {
            hash = ((hash << 5) + hash) + hash_array[j];
        }
    }
    
    // 处理剩余字节
    for (; i < len; i++) {
        hash = ((hash << 5) + hash) + str[i];
    }
    
    return hash % size;
}

// 缓存对齐的内存分配
static inline void* aligned_malloc(size_t size) {
    void* ptr = NULL;
    if (posix_memalign(&ptr, CACHE_LINE_SIZE, size)) {
        return NULL;
    }
    return ptr;
}

// 线程本地存储初始化
static inline void init_thread_local_storage(size_t size) {
    if (!thread_local_buffer || thread_local_buffer_size < size) {
        free(thread_local_buffer);
        thread_local_buffer = (char*)aligned_malloc(size);
        thread_local_buffer_size = size;
    }
}

// 优化的哈希表操作
HashMap* create_hashmap(int size) {
    HashMap* map = (HashMap*)aligned_malloc(sizeof(HashMap));
    map->size = size;
    map->buckets = (HashNode**)aligned_malloc(size * sizeof(HashNode*));
    memset(map->buckets, 0, size * sizeof(HashNode*));
    return map;
}

void hashmap_put(HashMap* map, const char* key, int position) {
    unsigned int index = optimized_hash(key, map->size);
    HashNode* node = map->buckets[index];
    
    // 查找已存在的键
    while (node != NULL) {
        if (avx2_strcmp(node->key, key, strlen(key))) {
            if (node->count >= node->capacity) {
                node->capacity *= 2;
                node->positions = (int*)aligned_malloc(sizeof(int) * node->capacity);
            }
            node->positions[node->count++] = position;
            return;
        }
        node = node->next;
    }
    
    // 创建新节点
    node = (HashNode*)aligned_malloc(sizeof(HashNode));
    node->key = strdup(key);
    node->capacity = 16;
    node->positions = (int*)aligned_malloc(sizeof(int) * node->capacity);
    node->positions[0] = position;
    node->count = 1;
    node->next = map->buckets[index];
    map->buckets[index] = node;
}

// 优化的重复查找函数
RepeatPattern* find_repeats(const char* query, const char* reference, int* repeat_count) {
    omp_set_num_threads(NUM_LOGICAL_CORES);
    
    int query_len = strlen(query);
    int ref_len = strlen(reference);
    
    printf("查询序列长度: %d\n", query_len);
    printf("参考序列长度: %d\n", ref_len);
    
    RepeatPattern* repeats = (RepeatPattern*)aligned_malloc(sizeof(RepeatPattern) * MAX_REPEATS);
    *repeat_count = 0;
    
    // 创建更大的哈希表以减少冲突
    HashMap* window_positions = create_hashmap(HASH_TABLE_SIZE);
    
    // 计算最大可能长度
    int max_possible_length = (query_len < MAX_LENGTH) ? query_len : MAX_LENGTH;
    
    // 按长度并行处理
    #pragma omp parallel
    {
        // 每个线程处理特定范围的长度
        #pragma omp for schedule(dynamic, 1)
        for (int length = MIN_LENGTH; length <= max_possible_length; length++) {
            if (length > query_len) continue;
            
            // 初始化线程本地存储
            init_thread_local_storage(length + 1);
            
            // 为当前长度构建查询序列的位置映射
            HashMap* local_positions = create_hashmap(HASH_TABLE_SIZE);
            
            // 分块处理查询序列
            for (int i = 0; i <= query_len - length; i += BLOCK_SIZE) {
                int end = i + BLOCK_SIZE;
                if (end > query_len - length + 1) end = query_len - length + 1;
                
                for (int j = i; j < end; j++) {
                    char* segment = thread_local_buffer;
                    memcpy(segment, query + j, length);
                    segment[length] = '\0';
                    
                    hashmap_put(local_positions, segment, j);
                }
            }
            
            // 并行处理参考序列
            for (int i = 0; i <= ref_len - length; i += BLOCK_SIZE) {
                int end = i + BLOCK_SIZE;
                if (end > ref_len - length + 1) end = ref_len - length + 1;
                
                for (int j = i; j < end; j++) {
                    
                    char* segment = thread_local_buffer;
                    memcpy(segment, reference + j, length);
                    segment[length] = '\0';
                    
                    // 检查正向重复
                    int pos_count;
                    int* positions = hashmap_get(local_positions, segment, &pos_count);
                    
                    if (positions && pos_count >= 2) {
                        process_matches(repeats, repeat_count, positions, pos_count, j, length, segment, false);
                    }
                    
                    // 检查反向互补重复
                    char* rev_comp = (char*)aligned_malloc(length + 1);
                    get_reverse_complement(segment, rev_comp, length);
                    
                    positions = hashmap_get(local_positions, rev_comp, &pos_count);
                    if (positions && pos_count >= 2) {
                        process_matches(repeats, repeat_count, positions, pos_count, j, length, segment, true);
                    }
                    
                    free(rev_comp);
                }
            }
            
            free_hashmap(local_positions);
        }
    }
    
    // 过滤和排序结果
    filter_nested_repeats(repeats, repeat_count);
    sort_repeats(repeats, *repeat_count);
    
    free_hashmap(window_positions);
    return repeats;
}

// AMD Ryzen 9 7940HX specific optimizations
#define NUM_PHYSICAL_CORES 16
#define NUM_LOGICAL_CORES 32
#define L1_CACHE_SIZE 32768     // 32KB per core
#define L2_CACHE_SIZE 512000    // 512KB per core
#define L3_CACHE_SIZE 32768000  // 32MB shared
#define CACHE_LINE_SIZE 64
#define MAX_LENGTH 120
#define MIN_LENGTH 10
#define MAX_REPEATS 1000

// Optimization macros
#define ALIGN_TO_CACHE __attribute__((aligned(CACHE_LINE_SIZE)))
#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

// 重复片段的数据结构
typedef struct ALIGN_TO_CACHE {
    int position;
    int length;
    int repeat_count;
    bool is_reverse;
    char* original_sequence;
    int query_position;
} RepeatPattern;

// 哈希表节点结构
typedef struct HashNode {
    char* key;
    int* positions;
    int count;
    int capacity;
    struct HashNode* next;
} HashNode ALIGN_TO_CACHE;

// 哈希表结构
typedef struct {
    HashNode** buckets;
    int size;
} HashMap ALIGN_TO_CACHE;

// Thread-local storage
static __thread char* thread_local_buffer = NULL;
static __thread size_t thread_local_buffer_size = 0;

// Function declarations
static inline void* aligned_malloc(size_t size) {
    void* ptr = NULL;
    if (posix_memalign(&ptr, CACHE_LINE_SIZE, size)) {
        return NULL;
    }
    return ptr;
}

static inline void init_thread_local_storage(size_t size) {
    if (!thread_local_buffer || thread_local_buffer_size < size) {
        free(thread_local_buffer);
        thread_local_buffer = (char*)aligned_malloc(size);
        thread_local_buffer_size = size;
    }
}

static inline void prefetch_read(const void* addr) {
    _mm_prefetch((const char*)addr, _MM_HINT_T0);
}

// Optimized string comparison using AVX2
static inline int avx2_strcmp(const char* str1, const char* str2, size_t len) {
    size_t i = 0;
    
    // Use AVX2 for blocks of 32 bytes
    for (; i + 32 <= len; i += 32) {
        __m256i v1 = _mm256_loadu_si256((__m256i*)(str1 + i));
        __m256i v2 = _mm256_loadu_si256((__m256i*)(str2 + i));
        __m256i vcmp = _mm256_cmpeq_epi8(v1, v2);
        int mask = _mm256_movemask_epi8(vcmp);
        if (mask != 0xFFFFFFFF) {
            return 0;
        }
    }
    
    // Handle remaining bytes
    return memcmp(str1 + i, str2 + i, len - i) == 0;
}

// 优化的哈希函数
static inline unsigned int optimized_hash(const char* str, int size) {
    unsigned int hash = 5381;
    int c;
    
    // Use SIMD for hash calculation when possible
    const __m256i mult = _mm256_set1_epi32(33);
    __m256i hash_vec = _mm256_set1_epi32(hash);
    
    // Process 32 bytes at a time
    while (strlen(str) >= 32) {
        __m256i data = _mm256_loadu_si256((__m256i*)str);
        hash_vec = _mm256_add_epi32(
            _mm256_mullo_epi32(hash_vec, mult),
            _mm256_cvtepu8_epi32(_mm256_extracti128_si256(data, 0))
        );
        str += 32;
    }
    
    // Combine hash values
    int* hash_array = (int*)&hash_vec;
    for (int i = 0; i < 8; i++) {
        hash = ((hash << 5) + hash) + hash_array[i];
    }
    
    // Process remaining bytes
    while ((c = *str++)) {
        hash = ((hash << 5) + hash) + c;
    }
    
    return hash % size;
}

// 函数声明
char* read_sequence(const char* filename);
void get_reverse_complement(const char* dna, char* result, int length);
int* build_next(const char* pattern, int length);
int* find_all_matches(const char* text, const char* pattern, int text_len, int pattern_len, int* count);
HashMap* create_hashmap(int size);
void hashmap_put(HashMap* map, const char* key, int position);
int* hashmap_get(HashMap* map, const char* key, int* count);
void clear_hashmap(HashMap* map);
void free_hashmap(HashMap* map);
RepeatPattern* find_repeats(const char* query, const char* reference, int* repeat_count);
void save_repeats_to_file(RepeatPattern* repeats, int count);
int* find_consecutive_groups(int* positions, int pos_count, int length, int* group_count);
void filter_nested_repeats(RepeatPattern* repeats, int* count);

// 添加缺失的函数声明
static void process_matches(RepeatPattern* repeats, int* repeat_count, 
                          int* positions, int pos_count, 
                          int position, int length, 
                          const char* sequence, bool is_reverse);
static void sort_repeats(RepeatPattern* repeats, int count);

// KMP算法实现
int* build_next(const char* pattern, int length) {
    int* next = (int*)malloc(sizeof(int) * length);
    next[0] = -1;
    int k = -1;
    int j = 0;
    
    while (j < length - 1) {
        if (k == -1 || pattern[j] == pattern[k]) {
            k++;
            j++;
            // 优化：检查下一个字符
            if (j < length && pattern[j] != pattern[k]) {
                next[j] = k;
            } else {
                next[j] = next[k];
            }
        } else {
            k = next[k];
        }
    }
    return next;
}

// 优化的查找所有匹配位置函数
int* find_all_matches(const char* text, const char* pattern, int text_len, int pattern_len, int* count) {
    int* next = (int*)aligned_malloc(sizeof(int) * pattern_len);
    int* matches = (int*)aligned_malloc(sizeof(int) * text_len);
    *count = 0;
    
    // Build KMP next array
    #pragma omp parallel num_threads(1)
    {
        next[0] = -1;
        int k = -1;
        int j = 0;
        
        while (j < pattern_len - 1) {
            if (k == -1 || pattern[j] == pattern[k]) {
                k++;
                j++;
                if (j < pattern_len && pattern[j] != pattern[k]) {
                    next[j] = k;
                } else {
                    next[j] = next[k];
                }
            } else {
                k = next[k];
            }
        }
    }
    
    // Parallel pattern matching
    #pragma omp parallel num_threads(NUM_LOGICAL_CORES)
    {
        int local_count = 0;
        int* local_matches = (int*)aligned_malloc(sizeof(int) * text_len);
        
        #pragma omp for schedule(guided)
        for (int start = 0; start < text_len; start += CACHE_LINE_SIZE) {
            int end = start + CACHE_LINE_SIZE;
            if (end > text_len) end = text_len;
            
            int i = start;
            int j = 0;
            
            while (i < end) {
                prefetch_read(text + i + CACHE_LINE_SIZE);
                
                if (j == -1 || avx2_strcmp(text + i, pattern + j, 1)) {
                    i++;
                    j++;
                } else {
                    j = next[j];
                }
                
                if (j == pattern_len) {
                    local_matches[local_count++] = i - j;
                    j = next[j-1];
                }
            }
        }
        
        // Merge local results
        #pragma omp critical
        {
            memcpy(matches + *count, local_matches, local_count * sizeof(int));
            *count += local_count;
        }
        
        free(local_matches);
    }
    
    free(next);
    return matches;
}

// 哈希表实现
HashMap* create_hashmap(int size) {
    HashMap* map = (HashMap*)malloc(sizeof(HashMap));
    map->size = size;
    map->buckets = (HashNode**)calloc(size, sizeof(HashNode*));
    return map;
}

// 计算哈希值
unsigned int hash_function(const char* str, int size) {
    unsigned int hash = 5381;
    int c;
    while ((c = *str++)) {
        hash = ((hash << 5) + hash) + c;
    }
    return hash % size;
}

void hashmap_put(HashMap* map, const char* key, int position) {
    unsigned int index = hash_function(key, map->size);
    HashNode* node = map->buckets[index];
    
    // 查找是否已存在该键
    while (node != NULL) {
        if (strcmp(node->key, key) == 0) {
            // 键存在，添加新位置
            if (node->count >= node->capacity) {
                node->capacity *= 2;
                node->positions = (int*)realloc(node->positions, sizeof(int) * node->capacity);
            }
            node->positions[node->count++] = position;
            return;
        }
        node = node->next;
    }
    
    // 创建新节点
    node = (HashNode*)malloc(sizeof(HashNode));
    node->key = strdup(key);
    node->capacity = 16;
    node->positions = (int*)malloc(sizeof(int) * node->capacity);
    node->positions[0] = position;
    node->count = 1;
    node->next = map->buckets[index];
    map->buckets[index] = node;
}

int* hashmap_get(HashMap* map, const char* key, int* count) {
    unsigned int index = hash_function(key, map->size);
    HashNode* node = map->buckets[index];
    
    while (node != NULL) {
        if (strcmp(node->key, key) == 0) {
            *count = node->count;
            return node->positions;
        }
        node = node->next;
    }
    
    *count = 0;
    return NULL;
}

void clear_hashmap(HashMap* map) {
    for (int i = 0; i < map->size; i++) {
        HashNode* node = map->buckets[i];
        while (node != NULL) {
            HashNode* next = node->next;
            free(node->key);
            free(node->positions);
            free(node);
            node = next;
        }
        map->buckets[i] = NULL;
    }
}

void free_hashmap(HashMap* map) {
    clear_hashmap(map);
    free(map->buckets);
    free(map);
}

// 读取序列
char* read_sequence(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("无法打开文件: %s\n", filename);
        exit(1);
    }
    
    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    fseek(file, 0, SEEK_SET);
    
    char* sequence = (char*)malloc(length + 1);
    size_t bytes_read = fread(sequence, 1, length, file);
    if (bytes_read != length) {
        printf("读取文件失败: %s\n", filename);
        free(sequence);
        fclose(file);
        exit(1);
    }
    sequence[length] = '\0';
    
    // 转换为大写并移除空白字符
    int write = 0;
    for (int read = 0; sequence[read]; read++) {
        char c = sequence[read];
        if (c >= 'a' && c <= 'z') {
            sequence[write++] = c - 32;
        } else if (c >= 'A' && c <= 'Z') {
            sequence[write++] = c;
        }
    }
    sequence[write] = '\0';
    
    fclose(file);
    return sequence;
}

// 获取反向互补序列
void get_reverse_complement(const char* dna, char* result, int length) {
    for (int i = 0; i < length; i++) {
        char c = dna[length - 1 - i];
        switch (c) {
            case 'A': result[i] = 'T'; break;
            case 'T': result[i] = 'A'; break;
            case 'G': result[i] = 'C'; break;
            case 'C': result[i] = 'G'; break;
            default: result[i] = 'N';
        }
    }
    result[length] = '\0';
}

// 查找连续重复组
int* find_consecutive_groups(int* positions, int pos_count, int length, int* group_count) {
    if (pos_count < 2) {
        *group_count = 0;
        return NULL;
    }
    
    int* groups = (int*)malloc(sizeof(int) * pos_count);
    int current_count = 1;
    *group_count = 0;
    
    for (int i = 1; i < pos_count; i++) {
        if (positions[i] == positions[i-1] + length) {
            current_count++;
        } else {
            if (current_count >= 2) {
                groups[(*group_count)++] = current_count;
            }
            current_count = 1;
        }
    }
    
    if (current_count >= 2) {
        groups[(*group_count)++] = current_count;
    }
    
    return groups;
}

// 过滤嵌套重复
void filter_nested_repeats(RepeatPattern* repeats, int* count) {
    if (*count <= 1) return;
    
    // 按位置和方向分组，保留最长的重复
    bool* to_remove = (bool*)calloc(*count, sizeof(bool));
    int new_count = *count;
    
    for (int i = 0; i < *count; i++) {
        if (to_remove[i]) continue;
        
        for (int j = i + 1; j < *count; j++) {
            if (to_remove[j]) continue;
            
            if (repeats[i].position == repeats[j].position && 
                repeats[i].is_reverse == repeats[j].is_reverse) {
                // 保留较长的重复
                if (repeats[i].length > repeats[j].length) {
                    to_remove[j] = true;
                    new_count--;
                } else {
                    to_remove[i] = true;
                    new_count--;
                    break;
                }
            }
        }
    }
    
    // 重新排列数组，移除被标记的重复
    int write = 0;
    for (int read = 0; read < *count; read++) {
        if (!to_remove[read]) {
            if (write != read) {
                repeats[write] = repeats[read];
            }
            write++;
        } else {
            free(repeats[read].original_sequence);
        }
    }
    
    *count = new_count;
    free(to_remove);
}

// 比较重复片段的函数
int compare_repeats(const void* a, const void* b) {
    RepeatPattern* repeatA = (RepeatPattern*)a;
    RepeatPattern* repeatB = (RepeatPattern*)b;
    int scoreA = repeatA->length * repeatA->repeat_count;
    int scoreB = repeatB->length * repeatB->repeat_count;
    return scoreB - scoreA; // 降序排列
}

// 优化的重复片段查找函数
RepeatPattern* find_repeats(const char* query, const char* reference, int* repeat_count) {
    // 初始化 OpenMP
    omp_set_num_threads(16);
    
    int query_len = strlen(query);
    int ref_len = strlen(reference);
    
    printf("查询序列长度: %d\n", query_len);
    printf("参考序列长度: %d\n", ref_len);
    
    RepeatPattern* repeats = (RepeatPattern*)aligned_malloc(sizeof(RepeatPattern) * MAX_REPEATS);
    *repeat_count = 0;
    
    // 创建哈希表用于存储片段位置
    HashMap* window_positions = create_hashmap(16384);
    
    // 特别关注区域
    int special_check_around = 400;
    
    // 按长度并行处理
    int max_possible_length = (query_len < MAX_LENGTH) ? query_len : MAX_LENGTH;
    #pragma omp parallel for schedule(dynamic) shared(repeats, repeat_count)
    for (int length = MIN_LENGTH; length <= max_possible_length; length++) {
        // 如果超出范围，跳过此次迭代
        if (length > query_len) continue;
        
        // 线程局部存储
        init_thread_local_storage(length + 1);
        
        HashMap* local_positions = create_hashmap(16384);
        
        // 构建查询序列的位置映射
        #pragma omp parallel for schedule(guided)
        for (int i = 0; i <= query_len - length; i++) {
            char* segment = thread_local_buffer;
            memcpy(segment, query + i, length);
            segment[length] = '\0';
            
            #pragma omp critical
            {
                hashmap_put(local_positions, segment, i);
            }
        }
        
        // 检查参考序列中的片段
        char* rev_comp = (char*)aligned_malloc(length + 1);
        
        #pragma omp parallel for schedule(guided)
        for (int i = 0; i <= ref_len - length; i++) {
            char* segment = thread_local_buffer;
            memcpy(segment, reference + i, length);
            segment[length] = '\0';
            
            // 使用AVX2优化的字符串比较
            int pos_count;
            int* positions = hashmap_get(local_positions, segment, &pos_count);
            
            if (positions && pos_count >= 2) {
                process_matches(repeats, repeat_count, positions, pos_count, i, length, segment, false);
            }
            
            // 检查反向互补重复
            get_reverse_complement(segment, rev_comp, length);
            positions = hashmap_get(local_positions, rev_comp, &pos_count);
            
            if (positions && pos_count >= 2) {
                process_matches(repeats, repeat_count, positions, pos_count, i, length, segment, true);
            }
        }
        
        free(rev_comp);
        free_hashmap(local_positions);
    }
    
    // 过滤和排序结果
    filter_nested_repeats(repeats, repeat_count);
    sort_repeats(repeats, *repeat_count);
    
    return repeats;
}

// 添加process_matches函数实现
static void process_matches(RepeatPattern* repeats, int* repeat_count, 
                          int* positions, int pos_count, 
                          int position, int length, 
                          const char* sequence, bool is_reverse) {
    int group_count;
    int* groups = find_consecutive_groups(positions, pos_count, length, &group_count);
    
    #pragma omp critical
    {
        for (int g = 0; g < group_count; g++) {
            if (*repeat_count >= MAX_REPEATS) break;
            
            repeats[*repeat_count].position = position+length;
            repeats[*repeat_count].length = length;
            repeats[*repeat_count].repeat_count = groups[g];
            repeats[*repeat_count].is_reverse = is_reverse;
            repeats[*repeat_count].original_sequence = strdup(sequence);
            repeats[*repeat_count].query_position = positions[g];
            (*repeat_count)++;
        }
    }
    
    free(groups);
}

// 添加sort_repeats函数实现
static void sort_repeats(RepeatPattern* repeats, int count) {
    // 使用快速排序，按长度和重复次数的乘积降序排序
    if (count < 2) return;

    RepeatPattern pivot = repeats[count / 2];
    int left = 0;
    int right = count - 1;

    while (left <= right) {
        while (repeats[left].length * repeats[left].repeat_count > 
               pivot.length * pivot.repeat_count) {
            left++;
        }
        while (repeats[right].length * repeats[right].repeat_count < 
               pivot.length * pivot.repeat_count) {
            right--;
        }
        if (left <= right) {
            RepeatPattern temp = repeats[left];
            repeats[left] = repeats[right];
            repeats[right] = temp;
            left++;
            right--;
        }
    }

    sort_repeats(repeats, right + 1);
    sort_repeats(repeats + left, count - left);
}

// 保存结果到文件
void save_repeats_to_file(RepeatPattern* repeats, int count) {
    FILE* file = fopen("repeat_results_new.txt", "w");
    if (!file) {
        printf("无法创建输出文件\n");
        return;
    }
    
    fprintf(file, "参考位置,长度,重复次数,是否反向重复,原始序列,查询位置\n");
    for (int i = 0; i < count; i++) {
        repeats[i].position += repeats[i].length;
        fprintf(file, "%d,%d,%d,%s,%s,%d\n",
                repeats[i].position,
                repeats[i].length,
                repeats[i].repeat_count,
                repeats[i].is_reverse ? "是" : "否",
                repeats[i].original_sequence,
                repeats[i].query_position);
    }
    fclose(file);
    
    // 保存详细信息
    file = fopen("repeat_details_new.txt", "w");
    if (!file) {
        printf("无法创建详细输出文件\n");
        return;
    }
    
    for (int i = 0; i < count; i++) {
        fprintf(file, "重复 #%d:\n", i+1);
        fprintf(file, "  参考位置: %d\n", repeats[i].position+repeats[i].length);
        fprintf(file, "  长度: %d\n", repeats[i].length);
        fprintf(file, "  重复次数: %d\n", repeats[i].repeat_count);
        fprintf(file, "  是否反向重复: %s\n", repeats[i].is_reverse ? "是" : "否");
        fprintf(file, "  原始序列: %s\n", repeats[i].original_sequence);
        fprintf(file, "  查询位置: %d\n\n", repeats[i].query_position);
    }
    fclose(file);
}

int main(int argc, char* argv[]) {
    // Initialize OpenMP
    omp_set_num_threads(NUM_LOGICAL_CORES);
    
    char* query_file = "query.txt";
    char* reference_file = "reference.txt";
    
    // Check command line arguments
    if (argc >= 3) {
        reference_file = argv[1];
        query_file = argv[2];
    }
    
    printf("读取查询序列: %s\n", query_file);
    char* query = read_sequence(query_file);
    
    printf("读取参考序列: %s\n", reference_file);
    char* reference = read_sequence(reference_file);
    
    // Start timing
    clock_t start = clock();
    
    // Find repeats
    int repeat_count;
    RepeatPattern* repeats = find_repeats(query, reference, &repeat_count);
    
    // Calculate execution time
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    //print results
    for (int i = 0; i < repeat_count; i++) {
        printf("position: %d, length: %d, repeat_count: %d, is_reverse: %d, original_sequence: %s, query_position: %d\n",
        repeats[i].position, repeats[i].length, repeats[i].repeat_count, repeats[i].is_reverse, repeats[i].original_sequence, repeats[i].query_position);
    }
    printf("找到 %d 个重复片段，耗时: %.2f ms\n", repeat_count, time_spent*1000);
    
    // Save results
    save_repeats_to_file(repeats, repeat_count);
    
    // Cleanup
    for (int i = 0; i < repeat_count; i++) {
        free(repeats[i].original_sequence);
    }
    free(repeats);
    free(query);
    free(reference);
    
    return 0;
}