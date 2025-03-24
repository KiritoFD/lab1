#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <memory>
#include <algorithm>

// 前向声明
class DNASequence;
class HashTable;
class RepeatFinder;
class FuzzyMatcher;

// DNA序列类
class DNASequence {
private:
    char* sequence;
    int length;

public:
    DNASequence(const char* filename) {
        sequence = readFromFile(filename);
        length = strlen(sequence);
    }

    ~DNASequence() {
        if (sequence) {
            free(sequence);
        }
    }

    const char* getSequence() const { return sequence; }
    int getLength() const { return length; }

    // 获取子序列
    char* getSubsequence(int start, int length) const {
        char* sub = (char*)malloc(length + 1);
        strncpy(sub, sequence + start, length);
        sub[length] = '\0';
        return sub;
    }

    // 获取反向互补序列
    static char* getReverseComplement(const char* seq, int length) {
        char* result = (char*)malloc(length + 1);
        for (int i = 0; i < length; i++) {
            char c = seq[length - 1 - i];
            switch (c) {
                case 'A': result[i] = 'T'; break;
                case 'T': result[i] = 'A'; break;
                case 'G': result[i] = 'C'; break;
                case 'C': result[i] = 'G'; break;
                default: result[i] = 'N';
            }
        }
        result[length] = '\0';
        return result;
    }

private:
    char* readFromFile(const char* filename) {
        FILE* file = fopen(filename, "r");
        if (!file) {
            printf("无法打开文件: %s\n", filename);
            exit(1);
        }

        fseek(file, 0, SEEK_END);
        long length = ftell(file);
        fseek(file, 0, SEEK_SET);

        char* seq = (char*)malloc(length + 1);
        size_t bytes_read = fread(seq, 1, length, file);
        if (bytes_read != length) {
            printf("读取文件失败: %s\n", filename);
            free(seq);
            fclose(file);
            exit(1);
        }
        seq[length] = '\0';

        // 转换为大写并移除空白字符
        int write = 0;
        for (int read = 0; seq[read]; read++) {
            char c = seq[read];
            if (c >= 'a' && c <= 'z') {
                seq[write++] = c - 32;
            } else if (c >= 'A' && c <= 'Z') {
                seq[write++] = c;
            }
        }
        seq[write] = '\0';

        fclose(file);
        return seq;
    }
};

// 模糊匹配类
class FuzzyMatcher {
private:
    float similarityThreshold;
    int windowSize;

public:
    FuzzyMatcher(float threshold = 0.85f, int wSize = 3) 
        : similarityThreshold(threshold), windowSize(wSize) {}

    bool isMatch(const char* str1, const char* str2, int length) const {
        return calculateSimilarity(str1, str2, length) >= similarityThreshold;
    }

    unsigned int getHash(const char* str, int length) const {
        unsigned int hash = 5381;
        for (int i = 0; i <= length - windowSize; i++) {
            unsigned int local_hash = 0;
            for (int j = 0; j < windowSize; j++) {
                local_hash = ((local_hash << 5) + local_hash) + str[i + j];
            }
            hash = ((hash << 5) + hash) + local_hash;
        }
        return hash;
    }

private:
    float calculateSimilarity(const char* seq1, const char* seq2, int length) const {
        int matches = 0;
        for (int i = 0; i < length; i++) {
            if (seq1[i] == seq2[i]) {
                matches++;
            }
        }
        return (float)matches / length;
    }
};

// 添加任务分配器类
class TaskDistributor {
private:
    int total_length;
    int min_segment_size;
    std::vector<std::pair<int, int>> segments;
    std::mutex mtx;
    size_t current_segment;

public:
    TaskDistributor(int total, int min_size = 1000) 
        : total_length(total), min_segment_size(min_size), current_segment(0) {
        // 根据序列长度动态计算段数
        int optimal_segments = std::thread::hardware_concurrency() * 2;
        int segment_size = total_length / optimal_segments;
        
        // 确保每段至少有最小长度
        if (segment_size < min_segment_size) {
            segment_size = min_segment_size;
            optimal_segments = (total_length + segment_size - 1) / segment_size;
        }

        // 创建段
        for (int i = 0; i < optimal_segments; ++i) {
            int start = i * segment_size;
            int end = (i == optimal_segments - 1) ? total_length : (i + 1) * segment_size;
            segments.emplace_back(start, end);
        }
    }

    bool getNextSegment(int& start, int& end) {
        std::lock_guard<std::mutex> lock(mtx);
        if (current_segment >= segments.size()) {
            return false;
        }
        start = segments[current_segment].first;
        end = segments[current_segment].second;
        ++current_segment;
        return true;
    }

    size_t getSegmentCount() const {
        return segments.size();
    }
};

// 修改 HashTable 类，使其线程安全
class HashTable {
private:
    struct Node {
        char* key;
        int* positions;
        int count;
        int capacity;
        Node* next;

        Node(const char* k, int pos) {
            key = strdup(k);
            capacity = 16;
            positions = (int*)malloc(sizeof(int) * capacity);
            positions[0] = pos;
            count = 1;
            next = nullptr;
        }

        ~Node() {
            free(key);
            free(positions);
        }
    };

    Node** buckets;
    int size;
    FuzzyMatcher* matcher;
    std::vector<std::unique_ptr<std::mutex>> bucket_mutexes; // 修改互斥锁存储方式

public:
    HashTable(int table_size = 16384, FuzzyMatcher* m = nullptr) 
        : size(table_size), matcher(m ? m : new FuzzyMatcher()) {
        buckets = (Node**)calloc(size, sizeof(Node*));
        // 初始化互斥锁
        bucket_mutexes.reserve(size);
        for (int i = 0; i < size; i++) {
            bucket_mutexes.push_back(std::make_unique<std::mutex>());
        }
    }

    ~HashTable() {
        clear();
        free(buckets);
        if (matcher) {
            delete matcher;
        }
    }

    void put(const char* key, int position) {
        unsigned int index = matcher->getHash(key, strlen(key)) % size;
        std::lock_guard<std::mutex> lock(*bucket_mutexes[index]);

        Node* node = buckets[index];
        while (node) {
            if (strlen(node->key) == strlen(key) && 
                matcher->isMatch(node->key, key, strlen(key))) {
                if (node->count >= node->capacity) {
                    node->capacity *= 2;
                    node->positions = (int*)realloc(node->positions, 
                                                  sizeof(int) * node->capacity);
                }
                node->positions[node->count++] = position;
                return;
            }
            node = node->next;
        }

        Node* new_node = new Node(key, position);
        new_node->next = buckets[index];
        buckets[index] = new_node;
    }

    int* get(const char* key, int* count) const {
        unsigned int index = matcher->getHash(key, strlen(key)) % size;
        Node* node = buckets[index];
        int total_count = 0;
        int* all_positions = nullptr;
        int capacity = 16;

        all_positions = (int*)malloc(sizeof(int) * capacity);

        while (node) {
            if (strlen(node->key) == strlen(key) &&
                matcher->isMatch(node->key, key, strlen(key))) {
                if (total_count + node->count > capacity) {
                    capacity = (total_count + node->count) * 2;
                    all_positions = (int*)realloc(all_positions, 
                                                sizeof(int) * capacity);
                }
                memcpy(all_positions + total_count, node->positions, 
                       sizeof(int) * node->count);
                total_count += node->count;
            }
            node = node->next;
        }

        if (total_count == 0) {
            free(all_positions);
            all_positions = nullptr;
        }

        *count = total_count;
        return all_positions;
    }

    void clear() {
        for (int i = 0; i < size; i++) {
            Node* node = buckets[i];
            while (node) {
                Node* next = node->next;
                delete node;
                node = next;
            }
            buckets[i] = nullptr;
        }
    }
};

// 重复模式结构
struct RepeatPattern {
    int position;
    int length;
    int repeat_count;
    bool is_reverse;
    char* original_sequence;
    int query_position;
};

// 修改 RepeatFinder 类
class RepeatFinder {
private:
    static const int MAX_REPEATS = 1000;
    static const int MIN_LENGTH = 10;
    static const int MAX_LENGTH = 120;
    static const int NUM_SEGMENTS = 4; // 将序列分成4段
    static const int BOUNDARY_SIZE = 120; // 边界重叠区域大小

    DNASequence* query;
    DNASequence* reference;
    HashTable* hashTable;
    FuzzyMatcher* matcher;
    std::mutex results_mutex;
    std::atomic<bool> should_terminate{false};

    struct SegmentInfo {
        int start_pos;
        int end_pos;
        int min_length;
        int max_length;
    };

public:
    RepeatFinder(const char* query_file, const char* reference_file) {
        query = new DNASequence(query_file);
        reference = new DNASequence(reference_file);
        matcher = new FuzzyMatcher();
        hashTable = new HashTable(16384, matcher);
    }

    ~RepeatFinder() {
        delete query;
        delete reference;
        delete hashTable;
        // matcher 由 hashTable 负责删除
    }

    RepeatPattern* findRepeats(int* repeat_count) {
        printf("Query sequence length: %d\n", query->getLength());
        printf("Reference sequence length: %d\n", reference->getLength());

        std::vector<std::unique_ptr<HashTable>> thread_hash_tables;
        std::vector<std::vector<RepeatPattern>> thread_results;
        std::vector<int> thread_counts;
        std::vector<std::thread> threads;

        // 创建任务分配器
        TaskDistributor length_distributor(MAX_LENGTH - MIN_LENGTH + 1);
        TaskDistributor pos_distributor(reference->getLength());

        // 初始化线程数据结构
        size_t num_threads = std::thread::hardware_concurrency();
        thread_results.resize(num_threads);
        thread_counts.resize(num_threads);
        
        // 启动工作线程
        for (size_t i = 0; i < num_threads; ++i) {
            thread_hash_tables.push_back(
                std::make_unique<HashTable>(16384, new FuzzyMatcher()));
            threads.emplace_back(&RepeatFinder::workerThread, this,
                               i, thread_hash_tables[i].get(),
                               std::ref(thread_results[i]),
                               std::ref(thread_counts[i]),
                               std::ref(length_distributor),
                               std::ref(pos_distributor));
        }

        // 等待所有线程完成
        for (auto& thread : threads) {
            thread.join();
        }

        // 合并结果
        RepeatPattern* repeats = (RepeatPattern*)malloc(
            sizeof(RepeatPattern) * MAX_REPEATS);
        *repeat_count = 0;

        for (size_t i = 0; i < num_threads; ++i) {
            const auto& thread_result = thread_results[i];
            int count = thread_counts[i];

            if (*repeat_count + count > MAX_REPEATS) {
                count = MAX_REPEATS - *repeat_count;
            }

            if (count > 0) {
                memcpy(repeats + *repeat_count, thread_result.data(),
                       count * sizeof(RepeatPattern));
                *repeat_count += count;
            }

            if (*repeat_count >= MAX_REPEATS) {
                break;
            }
        }

        // 过滤和排序结果
        filterDuplicateRepeats(repeats, repeat_count);
        filterNestedRepeats(repeats, repeat_count);
        quickSortRepeats(repeats, 0, *repeat_count - 1);

        return repeats;
    }

private:
    void workerThread(size_t thread_id, HashTable* local_hash_table,
                     std::vector<RepeatPattern>& local_results,
                     int& local_count,
                     TaskDistributor& length_distributor,
                     TaskDistributor& pos_distributor) {
        int length_start, length_end;
        while (length_distributor.getNextSegment(length_start, length_end)) {
            for (int length = length_start + MIN_LENGTH; 
                 length <= length_end + MIN_LENGTH && length <= query->getLength(); 
                 length++) {
                
                local_hash_table->clear();

                // 构建查询序列的哈希表
                for (int i = 0; i <= query->getLength() - length; i++) {
                    char* segment = query->getSubsequence(i, length);
                    local_hash_table->put(segment, i);
                    free(segment);
                }

                // 处理参考序列的不同段
                int pos_start, pos_end;
                while (pos_distributor.getNextSegment(pos_start, pos_end)) {
                    for (int i = pos_start; 
                         i <= pos_end - length && i <= reference->getLength() - length; 
                         i++) {
                        
                        if (should_terminate) return;

                        char* segment = reference->getSubsequence(i, length);

                        // 检查正向重复
                        int pos_count;
                        int* positions = local_hash_table->get(segment, &pos_count);

                        if (positions && pos_count >= 2) {
                            std::lock_guard<std::mutex> lock(results_mutex);
                            addRepeatPatternToVector(local_results, local_count,
                                                   i, length, positions, pos_count,
                                                   segment, false);
                            free(positions);
                        }

                        // 检查反向互补重复
                        char* rev_comp = DNASequence::getReverseComplement(segment, length);
                        positions = local_hash_table->get(rev_comp, &pos_count);

                        if (positions && pos_count >= 2) {
                            std::lock_guard<std::mutex> lock(results_mutex);
                            addRepeatPatternToVector(local_results, local_count,
                                                   i, length, positions, pos_count,
                                                   rev_comp, true);
                            free(positions);
                        }

                        free(segment);
                        free(rev_comp);

                        if (local_count >= MAX_REPEATS) {
                            should_terminate = true;
                            return;
                        }
                    }
                }
            }
        }
    }

    void addRepeatPatternToVector(std::vector<RepeatPattern>& results,
                                int& count,
                                int position, int length,
                                int* positions, int pos_count,
                                const char* sequence, bool is_reverse) {
        int group_count;
        int* groups = findConsecutiveGroups(positions, pos_count, length, &group_count);

        for (int g = 0; g < group_count; g++) {
            if (count >= MAX_REPEATS) {
                free(groups);
                return;
            }

            RepeatPattern pattern;
            pattern.position = position;
            pattern.length = length;
            pattern.repeat_count = groups[g];
            pattern.is_reverse = is_reverse;
            pattern.original_sequence = strdup(sequence);
            pattern.query_position = positions[g];
            
            results.push_back(pattern);
            count++;
        }

        free(groups);
    }

    void processSegment(RepeatPattern* repeats, int* repeat_count,
                       int min_length, int max_length,
                       int start_pos, int end_pos,
                       HashTable* localHashTable) {
        // 遍历当前长度范围
        for (int length = min_length; 
             length <= max_length && length <= query->getLength(); 
             length++) {
            localHashTable->clear();

            // 构建查询序列的哈希表
            for (int i = 0; i <= query->getLength() - length; i++) {
                char* segment = query->getSubsequence(i, length);
                localHashTable->put(segment, i);
                free(segment);
            }

            // 在参考序列的指定范围内查找重复
            for (int i = start_pos; i <= end_pos - length && i <= reference->getLength() - length; i++) {
                char* segment = reference->getSubsequence(i, length);

                // 检查正向重复
                int pos_count;
                int* positions = localHashTable->get(segment, &pos_count);

                if (positions && pos_count >= 2) {
                    addRepeatPattern(repeats, repeat_count, i, length,
                                   positions, pos_count, segment, false);
                    free(positions);
                }

                // 检查反向互补重复
                char* rev_comp = DNASequence::getReverseComplement(segment, length);
                positions = localHashTable->get(rev_comp, &pos_count);

                if (positions && pos_count >= 2) {
                    addRepeatPattern(repeats, repeat_count, i, length,
                                   positions, pos_count, rev_comp, true);
                    free(positions);
                }

                free(segment);
                free(rev_comp);
            }
        }
    }

    void filterDuplicateRepeats(RepeatPattern* repeats, int* count) {
        if (*count <= 1) return;

        std::vector<bool> to_remove(*count, false);
        int new_count = *count;

        // 使用并行处理比较重复
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < *count; i++) {
            if (to_remove[i]) continue;

            for (int j = i + 1; j < *count; j++) {
                if (to_remove[j]) continue;

                // 检查是否为重复的结果（考虑边界重叠区域）
                if (repeats[i].position == repeats[j].position &&
                    repeats[i].length == repeats[j].length &&
                    repeats[i].is_reverse == repeats[j].is_reverse) {
                    #pragma omp critical
                    {
                        to_remove[j] = true;
                        new_count--;
                    }
                }
            }
        }

        // 压缩结果数组
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
    }

    void addRepeatPattern(RepeatPattern* repeats, int* repeat_count,
                         int position, int length, int* positions, int pos_count,
                         const char* sequence, bool is_reverse) {
        int group_count;
        int* groups = findConsecutiveGroups(positions, pos_count, length, &group_count);

        for (int g = 0; g < group_count; g++) {
            if (*repeat_count >= MAX_REPEATS) break;

            repeats[*repeat_count].position = position;
            repeats[*repeat_count].length = length;
            repeats[*repeat_count].repeat_count = groups[g];
            repeats[*repeat_count].is_reverse = is_reverse;
            repeats[*repeat_count].original_sequence = strdup(sequence);
            repeats[*repeat_count].query_position = positions[g];
            (*repeat_count)++;
        }

        free(groups);
    }

    // 查找连续重复组
    int* findConsecutiveGroups(int* positions, int pos_count, 
                              int length, int* group_count) {
        if (pos_count < 2) {
            *group_count = 0;
            return nullptr;
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
    void filterNestedRepeats(RepeatPattern* repeats, int* count) {
        if (*count <= 1) return;

        bool* to_remove = (bool*)calloc(*count, sizeof(bool));
        int new_count = *count;

        for (int i = 0; i < *count; i++) {
            if (to_remove[i]) continue;

            for (int j = i + 1; j < *count; j++) {
                if (to_remove[j]) continue;

                if (repeats[i].position == repeats[j].position && 
                    repeats[i].is_reverse == repeats[j].is_reverse) {
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

    // 快速排序
        void quickSortRepeats(RepeatPattern* repeats, int left, int right) {
            if (left < right) {
                int pivot = partition(repeats, left, right);
                quickSortRepeats(repeats, left, pivot - 1);
                quickSortRepeats(repeats, pivot + 1, right);
            }
        }
    
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
    
    public:
        void saveResults(RepeatPattern* repeats, int repeat_count) {
            FILE* output = fopen("repeat_results.txt", "w");
            if (!output) {
                printf("Error: Could not open output file for writing.\n");
                return;
            }
    
            fprintf(output, "Found %d repeat fragments\n", repeat_count);
            for (int i = 0; i < repeat_count; i++) {
                fprintf(output, "Repeat #%d: Position %d, Length %d, Repeat Count %d, Is Reverse Repeat %s\n",
                    i + 1,
                    repeats[i].position,
                    repeats[i].length,
                    repeats[i].repeat_count,
                    repeats[i].is_reverse ? "Yes" : "No");
            }
    
            fclose(output);
            printf("Results saved to repeat_results.txt\n");
        }
};

int main(int argc, char* argv[]) {
    const char* query_file = "query.txt"; // 修改为 const char*
    const char* reference_file = "reference.txt";

    if (argc >= 3) {
        reference_file = argv[1];
        query_file = argv[2];
    }

    // 开始计时
    clock_t start_time = clock();

    // 创建重复查找器并处理序列
    RepeatFinder finder(query_file, reference_file);
    int repeat_count;
    RepeatPattern* repeats = finder.findRepeats(&repeat_count);

    // 计算耗时
    clock_t end_time = clock();
    double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // 输出结果
    printf("Found %d repeat fragments, elapsed time: %.2f ms\n", 
           repeat_count, time_spent * 1000);
    for (int i = 0; i < repeat_count; i++) {
        printf("Repeat #%d: Position %d, Length %d, Repeat Count %d, "
               "Is Reverse Repeat %s\n",
               i + 1,
               repeats[i].position,
               repeats[i].length,
               repeats[i].repeat_count,
               repeats[i].is_reverse ? "Yes" : "No");
    }

    // 保存结果
    finder.saveResults(repeats, repeat_count);

    // 释放内存
    for (int i = 0; i < repeat_count; i++) {
        free(repeats[i].original_sequence);
    }
    free(repeats);

    return 0;
}