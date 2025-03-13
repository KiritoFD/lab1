#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#define MAX_LENGTH 120
#define MIN_LENGTH 10
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

// 查找重复片段主函数
RepeatPattern* find_repeats(const char* query, const char* reference, int* repeat_count) {
    int query_len = strlen(query);
    int ref_len = strlen(reference);
    
    printf("查询序列长度: %d\n", query_len);
    printf("参考序列长度: %d\n", ref_len);
    
    RepeatPattern* repeats = (RepeatPattern*)malloc(sizeof(RepeatPattern) * MAX_REPEATS);
    *repeat_count = 0;
    
    // 创建哈希表用于存储片段位置
    HashMap* window_positions = create_hashmap(16384); // 使用较大的哈希表以减少冲突
    
    // 特别关注区域
    int special_check_around = 400;
    
    // 按长度遍历
    for (int length = MIN_LENGTH; length <= MAX_LENGTH && length <= query_len; length++) {
        //printf("处理长度: %d\n", length);
        clear_hashmap(window_positions);
        
        // 为当前长度构建查询序列的位置映射
        char* segment = (char*)malloc(length + 1);
        for (int i = 0; i <= query_len - length; i++) {
            strncpy(segment, query + i, length);
            segment[length] = '\0';
            hashmap_put(window_positions, segment, i);
        }
        
        // 检查参考序列中的片段
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
    }
    
    // 过滤嵌套重复并排序
    filter_nested_repeats(repeats, repeat_count);
    
    // 使用快速排序优化重复片段排序（按长度和重复次数）
    qsort(repeats, *repeat_count, sizeof(RepeatPattern), compare_repeats);
    
    free_hashmap(window_positions);
    return repeats;
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
    char* query_file = "query.txt";
    char* reference_file = "reference.txt";
    
    // 检查命令行参数
    if (argc >= 3) {
        reference_file = argv[1];
        query_file = argv[2];
    }
    
    printf("读取查询序列: %s\n", query_file);
    char* query = read_sequence(query_file);
    
    printf("读取参考序列: %s\n", reference_file);
    char* reference = read_sequence(reference_file);
    
    // 开始计时
    clock_t start_time = clock();
    
    // 查找重复
    int repeat_count;
    RepeatPattern* repeats = find_repeats(query, reference, &repeat_count);
    
    // 计算耗时
    clock_t end_time = clock();
    double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    
    printf("找到 %d 个重复片段，耗时: %.2f 秒\n", repeat_count, time_spent);
    
    //打印结果
    for (int i = 0; i < repeat_count; i++) {
        printf("重复 #%d: 位置: %d, 长度: %d, 重复次数: %d, 是否反向: %s\n",
               i+1, repeats[i].position+repeats[i].length, repeats[i].length, repeats[i].repeat_count,
               repeats[i].is_reverse ? "是" : "否");
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