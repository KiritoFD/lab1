#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
//编译命令：gcc -march=znver4 -mtune=znver4 -O3 -ffast-math -flto -fuse-linker-plugin     -fprefetch-loop-arrays -funroll-loops -fomit-frame-pointer -mavx2 -mfma     -pthread -fopenmp     -DCPU_RYZEN_7940HX -DNUM_CORES=16 -DNUM_THREADS=32     -DL1_CACHE_SIZE=32768 -DL2_CACHE_SIZE=512000 -DL3_CACHE_SIZE=32768000     dna_repeat_finder_new.c -o dna_repeat_finder_new

//优化版编译命令：gcc -march=znver4 -mtune=znver4 -Ofast -flto -fuse-linker-plugin -fgraphite-identity -floop-nest-optimize -fprefetch-loop-arrays -funroll-loops -funroll-all-loops -fomit-frame-pointer -mavx2 -mfma  -pthread -fopenmp -fopt-info-vec-optimized -fmodulo-sched -fmodulo-sched-allow-regmoves -floop-interchange -floop-unroll-and-jam -ftree-loop-distribution -ftree-vectorize -funsafe-math-optimizations -ftracer -fweb -frename-registers -finline-functions -fipa-pta -falign-functions=64 -DCPU_RYZEN_7940HX -DNUM_CORES=16 -DNUM_THREADS=32 -DL1_CACHE_SIZE=32768 -DL2_CACHE_SIZE=512000 -DL3_CACHE_SIZE=32768000 dna_repeat_finder_new.c -o dna_repeat_finder_new

#define MAX_LENGTH 101
#define MIN_LENGTH 50
#define MAX_REPEATS 1000

// 重复片段的数据结构
typedef struct {
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
} HashNode;

// 哈希表结构
typedef struct {
    HashNode** buckets;
    int size;
    float similarity_threshold; // 相似度阈值
} HashMap;

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
void quick_sort_repeats(RepeatPattern* repeats, int left, int right);
HashMap* build_sequence_hashmap(const char* sequence, int seq_len, int length);

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

// 查找所有匹配位置
int* find_all_matches(const char* text, const char* pattern, int text_len, int pattern_len, int* count) {
    int* next = build_next(pattern, pattern_len);
    int* matches = (int*)malloc(sizeof(int) * text_len); // 最多text_len个匹配
    *count = 0;
    
    int i = 0; // 文本指针
    int j = 0; // 模式串指针
    
    while (i < text_len) {
        if (j == -1 || text[i] == pattern[j]) {
            i++;
            j++;
        } else {
            j = next[j];
        }
        
        if (j == pattern_len) {
            matches[*count] = i - j;
            (*count)++;
            j = next[j-1];
        }
    }
    
    free(next);
    return matches;
}

// 哈希表实现
HashMap* create_hashmap(int size) {
    HashMap* map = (HashMap*)malloc(sizeof(HashMap));
    map->size = size;
    map->buckets = (HashNode**)calloc(size, sizeof(HashNode*));
    map->similarity_threshold = 0.85f; // 默认相似度阈值
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

// DNA序列相似度计算
float calculate_similarity(const char* seq1, const char* seq2, int length) {
    int matches = 0;
    for (int i = 0; i < length; i++) {
        if (seq1[i] == seq2[i]) {
            matches++;
        }
    }
    return (float)matches / length;
}

// 获取模糊哈希值 - 使用滚动哈希提高效率
unsigned int fuzzy_hash_function(const char* str, int length, int size) {
    const unsigned int prime = 31; // 用于滚动哈希的素数
    unsigned int hash = 0;
    unsigned int power = 1; // prime^(window_size-1)
    int window_size = 3; // k-mer大小
    
    if (length < window_size) {
        // 如果字符串长度小于窗口大小，使用普通哈希
        for (int i = 0; i < length; i++) {
            hash = prime * hash + str[i];
        }
        return hash % size;
    }
    
    // 计算第一个窗口的哈希值
    for (int i = 0; i < window_size; i++) {
        hash = prime * hash + str[i];
    }
    
    // 计算prime^(window_size-1)，用于之后移除最左边字符的贡献
    for (int i = 0; i < window_size - 1; i++) {
        power = (power * prime);
    }
    
    unsigned int global_hash = hash;
    
    // 使用滚动哈希计算其余窗口的哈希值
    for (int i = window_size; i < length; i++) {
        // 移除最左边字符的贡献
        hash = hash - power * str[i - window_size];
        // 将哈希值乘以prime并添加新字符
        hash = hash * prime + str[i];
        // 将当前窗口的哈希值累加到全局哈希值
        global_hash = global_hash ^ hash; // 使用异或操作合并哈希值
    }
    
    return global_hash % size;
}

void hashmap_put(HashMap* map, const char* key, int position) {
    unsigned int index = fuzzy_hash_function(key, strlen(key), map->size);
    HashNode* node = map->buckets[index];
    
    // 查找是否存在相似的键
    while (node != NULL) {
        if (strlen(node->key) == strlen(key) && 
            calculate_similarity(node->key, key, strlen(key)) >= map->similarity_threshold) {
            // 找到相似的键，添加新位置
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
    unsigned int index = fuzzy_hash_function(key, strlen(key), map->size);
    HashNode* node = map->buckets[index];
    int total_count = 0;
    int* all_positions = NULL;
    int capacity = 16;
    
    all_positions = (int*)malloc(sizeof(int) * capacity);
    
    while (node != NULL) {
        if (strlen(node->key) == strlen(key) &&
            calculate_similarity(node->key, key, strlen(key)) >= map->similarity_threshold) {
            // 找到相似的键，合并位置
            if (total_count + node->count > capacity) {
                capacity = (total_count + node->count) * 2;
                all_positions = (int*)realloc(all_positions, sizeof(int) * capacity);
            }
            memcpy(all_positions + total_count, node->positions, sizeof(int) * node->count);
            total_count += node->count;
        }
        node = node->next;
    }
    
    if (total_count == 0) {
        free(all_positions);
        all_positions = NULL;
    }
    
    *count = total_count;
    return all_positions;
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

// 快速排序分区函数
int partition(RepeatPattern* repeats, int left, int right) {
    int pivot_value = repeats[right].length * repeats[right].repeat_count;
    int i = left - 1;
    
    for (int j = left; j < right; j++) {
        if (repeats[j].length * repeats[j].repeat_count >= pivot_value) {
            i++;
            RepeatPattern temp = repeats[i];
            repeats[i] = repeats[j];
            repeats[j] = temp;
        }
    }
    
    RepeatPattern temp = repeats[i + 1];
    repeats[i + 1] = repeats[right];
    repeats[right] = temp;
    
    return i + 1;
}

// 快速排序函数
void quick_sort_repeats(RepeatPattern* repeats, int left, int right) {
    if (left < right) {
        int pivot = partition(repeats, left, right);
        
        quick_sort_repeats(repeats, left, pivot - 1);
        quick_sort_repeats(repeats, pivot + 1, right);
    }
}

// 查找重复片段主函数
RepeatPattern* find_repeats(const char* query, const char* reference, int* repeat_count) {
    int query_len = strlen(query);
    int ref_len = strlen(reference);
    
    printf("Query sequence length: %d\n", query_len);
    printf("Reference sequence length: %d\n", ref_len);
    
    RepeatPattern* repeats = (RepeatPattern*)malloc(sizeof(RepeatPattern) * MAX_REPEATS);
    *repeat_count = 0;
    
    // 特别关注区域
    // 特别关注区域位置
    int special_check_around = query_len/2; // 默认值
    
    // 从环境变量获取值 (如果存在)
    char* env_value = getenv("SPECIAL_CHECK_AREA");
    if (env_value != NULL) {
        int parsed_value = atoi(env_value);
        if (parsed_value > 0) {
            special_check_around = parsed_value;
        }
    }
    
    printf("Special check area around position: %d\n", special_check_around);
    
    // 按长度遍历
    for (int length = MIN_LENGTH; length <= MAX_LENGTH && length <= query_len; length++) {
        // 为当前长度构建查询序列的位置映射
        HashMap* window_positions = build_sequence_hashmap(query, query_len, length);
        
        // 检查参考序列中的片段
        char* segment = (char*)malloc(length + 1);
        char* rev_comp = (char*)malloc(length + 1);
        
        for (int i = 0; i <= ref_len - length; i++) {
            // 对于非特殊区域限制检查长度范围
            if (abs(i - special_check_around) > 10 && length > MIN_LENGTH + 10) {
                if (length > 100) break;
            }
            
            strncpy(segment, reference + i, length);
            segment[length] = '\0';
            
            // 检查正向重复
            int pos_count;
            int* positions = hashmap_get(window_positions, segment, &pos_count);
            
            if (positions != NULL && pos_count >= 2) {
                int group_count;
                int* groups = find_consecutive_groups(positions, pos_count, length, &group_count);
                
                for (int g = 0; g < group_count; g++) {
                    if (*repeat_count >= MAX_REPEATS) break;
                    
                    repeats[*repeat_count].position = i;
                    repeats[*repeat_count].length = length;
                    repeats[*repeat_count].repeat_count = groups[g];
                    repeats[*repeat_count].is_reverse = false;
                    repeats[*repeat_count].original_sequence = strdup(segment);
                    repeats[*repeat_count].query_position = positions[g];
                    (*repeat_count)++;
                }
                
                free(groups);
            }
            
            // 检查反向互补重复
            get_reverse_complement(segment, rev_comp, length);
            positions = hashmap_get(window_positions, rev_comp, &pos_count);
            
            if (positions != NULL && pos_count >= 2) {
                int group_count;
                int* groups = find_consecutive_groups(positions, pos_count, length, &group_count);
                
                for (int g = 0; g < group_count; g++) {
                    if (*repeat_count >= MAX_REPEATS) break;
                    
                    repeats[*repeat_count].position = i;
                    repeats[*repeat_count].length = length;
                    repeats[*repeat_count].repeat_count = groups[g];
                    repeats[*repeat_count].is_reverse = true;
                    repeats[*repeat_count].original_sequence = strdup(segment);
                    repeats[*repeat_count].query_position = positions[g];
                    (*repeat_count)++;
                }
                
                free(groups);
            }
        }
        
        free(segment);
        free(rev_comp);
        free_hashmap(window_positions);
    }
    
    // 过滤嵌套重复并排序
    filter_nested_repeats(repeats, repeat_count);
    quick_sort_repeats(repeats, 0, *repeat_count - 1);
    
    return repeats;
}

// 保存结果到文件
void save_repeats_to_file(RepeatPattern* repeats, int count) {
    FILE* file = fopen("repeat_results_new.txt", "w");
    if (!file) {
        printf("Unable to create output file\n");
        return;
    }
    
    fprintf(file, "Reference Position,Length,Repeat Count,Is Reverse Repeat,Original Sequence,Query Position\n");
    for (int i = 0; i < count; i++) {
        repeats[i].position += repeats[i].length;
        fprintf(file, "%d,%d,%d,%s,%s,%d\n",
                repeats[i].position,
                repeats[i].length,
                repeats[i].repeat_count,
                repeats[i].is_reverse ? "Yes" : "No",
                repeats[i].original_sequence,
                repeats[i].query_position);
    }
    fclose(file);
    
    // Save detailed information
    file = fopen("repeat_details_new.txt", "w");
    if (!file) {
        printf("Unable to create detailed output file\n");
        return;
    }
    
    for (int i = 0; i < count; i++) {
        fprintf(file, "Repeat #%d:\n", i+1);
        fprintf(file, "  Reference Position: %d\n", repeats[i].position);
        fprintf(file, "  Length: %d\n", repeats[i].length);
        fprintf(file, "  Repeat Count: %d\n", repeats[i].repeat_count);
        fprintf(file, "  Is Reverse Repeat: %s\n", repeats[i].is_reverse ? "Yes" : "No");
        fprintf(file, "  Original Sequence: %s\n", repeats[i].original_sequence);
        fprintf(file, "  Query Position: %d\n\n", repeats[i].query_position);
    }
    fclose(file);
}

// 建立序列哈希表
HashMap* build_sequence_hashmap(const char* sequence, int seq_len, int length) {
    HashMap* window_positions = create_hashmap(16384); // 使用较大的哈希表以减少冲突
    
    // 构建查询序列的位置映射
    char* segment = (char*)malloc(length + 1);
    for (int i = 0; i <= seq_len - length; i++) {
        strncpy(segment, sequence + i, length);
        segment[length] = '\0';
        hashmap_put(window_positions, segment, i);
    }
    free(segment);
    
    return window_positions;
}

int main(int argc, char* argv[]) {
    char* query_file = "query.txt";
    char* reference_file = "reference.txt";
    
    // Check command-line arguments
    if (argc >= 3) {
        reference_file = argv[1];
        query_file = argv[2];
    }
    
    printf("Reading query sequence: %s\n", query_file);
    char* query = read_sequence(query_file);
    
    printf("Reading reference sequence: %s\n", reference_file);
    char* reference = read_sequence(reference_file);
    
    // Start timing
    clock_t start_time = clock();
    
    // Find repeats
    int repeat_count;
    RepeatPattern* repeats = find_repeats(query, reference, &repeat_count);
    
    // Calculate elapsed time
    clock_t end_time = clock();
    double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    
    printf("Found %d repeat fragments, elapsed time: %.2f ms\n", repeat_count, time_spent * 1000);
    // Print results
    for (int i = 0; i < repeat_count; i++) {
        printf("Repeat #%d: Position %d, Length %d, Repeat Count %d, Is Reverse Repeat %s\n",
               i + 1,
               repeats[i].position,
               repeats[i].length,
               repeats[i].repeat_count,
               repeats[i].is_reverse ? "Yes" : "No");
    }
    // 保存结果
    save_repeats_to_file(repeats, repeat_count);
    
    // 释放内存
    for (int i = 0; i < repeat_count; i++) {
        free(repeats[i].original_sequence);
    }
    free(repeats);
    free(query);
    free(reference);
    
    return 0;
}