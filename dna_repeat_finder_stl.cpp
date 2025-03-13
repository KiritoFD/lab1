/*
针对 AMD Ryzen 9 7940HX (Zen 4) 的优化编译命令:
g++ -std=c++20 -march=znver4 -mtune=znver4 -O3 -ffast-math -flto \
    -fuse-linker-plugin -fprefetch-loop-arrays -funroll-loops \
    -fomit-frame-pointer -mavx512f -mavx512dq -mavx512vl -mavx512bw \
    -pthread -fopenmp -ftree-vectorize -fopt-info-vec \
    dna_repeat_finder_stl.cpp -o dna_repeat_finder_stl
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <thread>
#include <mutex>
#include <functional>
#include <cstring>
#include <immintrin.h>
#include <omp.h>
#include <atomic>
#include <condition_variable>
#include <queue>

// AMD Ryzen 9 7940HX 特定参数
constexpr int L1_CACHE_SIZE = 32 * 1024;  // 32KB per core
constexpr int L2_CACHE_SIZE = 1024 * 1024;  // 1MB per core
constexpr int L3_CACHE_SIZE = 32 * 1024 * 1024;  // 32MB shared
constexpr int CACHE_LINE_SIZE = 64;
constexpr int NUM_PHYSICAL_CORES = 16;
constexpr int NUM_LOGICAL_CORES = 32;
constexpr int PREFETCH_DISTANCE = 256; // 预取距离优化

// 算法参数
constexpr int MIN_LENGTH = 10;
constexpr int MAX_LENGTH = 120;
constexpr size_t BLOCK_SIZE = 256; // 优化块大小以适应L1缓存

// 确保数据结构按缓存线对齐
alignas(CACHE_LINE_SIZE) alignas(32) static const char complement_table[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 'T', 0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0,
    0, 0, 0, 0, 'A', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 't', 0, 'g', 0, 0, 0, 'c', 0, 0, 0, 0, 0, 0, 'n', 0,
    0, 0, 0, 0, 'a', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// 重复片段的数据结构
struct alignas(CACHE_LINE_SIZE) RepeatPattern {
    int position;
    int length;
    int repeat_count;
    bool is_reverse;
    std::string original_sequence;
    int query_position;
    
    bool operator<(const RepeatPattern& other) const {
        return length * repeat_count > other.length * other.repeat_count;
    }
};

// 使用AVX2/AVX-512指令优化的字符串比较
inline bool simd_strcmp(const char* str1, const char* str2, size_t len) {
    size_t i = 0;
    
    // 预取数据到缓存
    _mm_prefetch(str1 + PREFETCH_DISTANCE, _MM_HINT_T0);
    _mm_prefetch(str2 + PREFETCH_DISTANCE, _MM_HINT_T0);
    
#ifdef __AVX512F__
    // 使用AVX-512处理64字节块
    for (; i + 64 <= len; i += 64) {
        __m512i v1 = _mm512_loadu_si512((const void*)(str1 + i));
        __m512i v2 = _mm512_loadu_si512((const void*)(str2 + i));
        __mmask64 mask = _mm512_cmpeq_epi8_mask(v1, v2);
        if (mask != 0xFFFFFFFFFFFFFFFFULL) {
            return false;
        }
        _mm_prefetch(str1 + i + PREFETCH_DISTANCE, _MM_HINT_T0);
        _mm_prefetch(str2 + i + PREFETCH_DISTANCE, _MM_HINT_T0);
    }
#endif

    // 使用AVX2处理32字节块
    for (; i + 32 <= len; i += 32) {
        __m256i v1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(str1 + i));
        __m256i v2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(str2 + i));
        __m256i vcmp = _mm256_cmpeq_epi8(v1, v2);
        int mask = _mm256_movemask_epi8(vcmp);
        if (mask != 0xFFFFFFFF) {
            return false;
        }
    }
    
    // 处理剩余字节
    return memcmp(str1 + i, str2 + i, len - i) == 0;
}

// 优化的反向互补序列生成
std::string get_reverse_complement(const std::string& dna) {
    const size_t len = dna.length();
    std::string result(len, 'N');
    
    // 并行处理反向互补
    #pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < len; ++i) {
        result[i] = complement_table[static_cast<unsigned char>(dna[len - 1 - i])];
    }
    
    return result;
}

// 任务定义
struct Task {
    int length;
    int start_pos;
    int end_pos;
    bool is_reference; // true=处理参考序列, false=处理查询序列
};

// 线程安全任务队列
class TaskQueue {
private:
    std::queue<Task> tasks;
    std::mutex mutex;
    std::condition_variable cv;
    bool done = false;

public:
    bool pop(Task& task) {
        std::unique_lock<std::mutex> lock(mutex);
        cv.wait(lock, [this] { return !tasks.empty() || done; });
        
        if (tasks.empty()) return false;
        
        task = std::move(tasks.front());
        tasks.pop();
        return true;
    }

    void push(Task&& task) {
        {
            std::unique_lock<std::mutex> lock(mutex);
            tasks.push(std::move(task));
        }
        cv.notify_one();
    }

    void finish() {
        {
            std::unique_lock<std::mutex> lock(mutex);
            done = true;
        }
        cv.notify_all();
    }
    
    size_t size() {
        std::unique_lock<std::mutex> lock(mutex);
        return tasks.size();
    }
};

// 全局变量定义
std::mutex g_io_mutex; // 用于输出的互斥锁
const std::string* g_query_ptr = nullptr;
const std::string* g_reference_ptr = nullptr;
std::mutex g_seq_mutex;

// 获取全局序列引用的函数
const std::string& get_query() {
    if (!g_query_ptr) {
        throw std::runtime_error("Query sequence not set");
    }
    return *g_query_ptr;
}

const std::string& get_reference() {
    if (!g_reference_ptr) {
        throw std::runtime_error("Reference sequence not set");
    }
    return *g_reference_ptr;
}

void set_sequences(const std::string& query, const std::string& reference) {
    std::lock_guard<std::mutex> lock(g_seq_mutex);
    g_query_ptr = &query;
    g_reference_ptr = &reference;
}

// 任务处理类
class TaskProcessor {
private:
    std::mutex positions_mutex;
    std::mutex results_mutex;
    std::unordered_map<std::string, std::vector<int>> positions;
    std::vector<RepeatPattern> results;
    
public:
    void process_query_segment(int length, int start_pos, int end_pos) {
        try {
            const std::string& query = get_query();
            const int query_len = query.length();
            
            // 使用线程本地存储减少锁竞争
            std::unordered_map<std::string, std::vector<int>> local_positions;
            local_positions.reserve((end_pos - start_pos) / length);
            
            for (int i = start_pos; i < end_pos && i <= query_len - length; ++i) {
                if ((i & 15) == 0) {
                    _mm_prefetch(query.data() + i + PREFETCH_DISTANCE, _MM_HINT_T0);
                }
                
                std::string segment = query.substr(i, length);
                local_positions[segment].push_back(i);
            }
            
            // 批量更新全局positions
            {
                std::lock_guard<std::mutex> lock(positions_mutex);
                for (auto& [segment, pos_vec] : local_positions) {
                    auto& target = positions[segment];
                    target.insert(target.end(), pos_vec.begin(), pos_vec.end());
                }
            }
        } catch (const std::exception& e) {
            std::lock_guard<std::mutex> lock(g_io_mutex);
            std::cerr << "处理查询序列段错误: " << e.what() << std::endl;
        }
    }
    
    void process_reference_segment(int length, int start_pos, int end_pos) {
        try {
            const std::string& reference = get_reference();
            const int ref_len = reference.length();
            std::vector<RepeatPattern> local_results;
            local_results.reserve(100);
            
            for (int i = start_pos; i < end_pos && i <= ref_len - length; ++i) {
                if (i % 5000 == 0) {
                    std::lock_guard<std::mutex> lock(g_io_mutex);
                    std::cout << "处理长度 " << length << " 进度: " 
                             << (float)i/ref_len*100.0f << "%\r" << std::flush;
                }
                
                if ((i & 15) == 0) {
                    _mm_prefetch(reference.data() + i + PREFETCH_DISTANCE, _MM_HINT_T0);
                }
                
                std::string segment = reference.substr(i, length);
                check_repeats(segment, i, length, false, local_results);
                
                std::string rev_comp = get_reverse_complement(segment);
                check_repeats(rev_comp, i, length, true, local_results);
            }
            
            if (!local_results.empty()) {
                std::lock_guard<std::mutex> lock(results_mutex);
                results.insert(results.end(), 
                             std::make_move_iterator(local_results.begin()),
                             std::make_move_iterator(local_results.end()));
            }
        } catch (const std::exception& e) {
            std::lock_guard<std::mutex> lock(g_io_mutex);
            std::cerr << "处理参考序列段错误: " << e.what() << std::endl;
        }
    }
    
    void check_repeats(const std::string& segment, int pos, int length, bool is_reverse,
                      std::vector<RepeatPattern>& local_results) {
        std::vector<int> pos_vec;
        {
            std::lock_guard<std::mutex> lock(positions_mutex);
            auto it = positions.find(segment);
            if (it == positions.end() || it->second.size() < 2) {
                return;
            }
            pos_vec = it->second; // 复制到本地处理
        }
        
        std::vector<std::vector<int>> consecutive_groups;
        if (!pos_vec.empty()) {
            std::vector<int> current_group = {pos_vec[0]};
            
            for (size_t k = 1; k < pos_vec.size(); ++k) {
                if (pos_vec[k] == current_group.back() + length) {
                    current_group.push_back(pos_vec[k]);
                } else {
                    if (current_group.size() >= 2) {
                        consecutive_groups.push_back(std::move(current_group));
                        current_group = {pos_vec[k]};
                    } else {
                        current_group[0] = pos_vec[k];
                    }
                }
            }
            
            if (current_group.size() >= 2) {
                consecutive_groups.push_back(std::move(current_group));
            }
        }
        
        for (const auto& group : consecutive_groups) {
            local_results.push_back({
                pos,
                length,
                static_cast<int>(group.size()),
                is_reverse,
                segment,
                group[0]
            });
        }
    }
    
    void clear_positions() {
        std::lock_guard<std::mutex> lock(positions_mutex);
        positions.clear();
    }
    
    std::vector<RepeatPattern>& get_results() {
        return results;
    }
};

// 优化的查找重复片段函数
std::vector<RepeatPattern> find_repeats(const std::string& query, const std::string& reference) {
    const int query_len = query.length();
    const int ref_len = reference.length();
    
    std::cout << "查询序列长度: " << query_len << std::endl;
    std::cout << "参考序列长度: " << ref_len << std::endl;
    
    // 设置全局序列引用
    set_sequences(query, reference);
    
    // 使用较少的线程数以减少竞争
    int optimal_threads = std::max(1, NUM_LOGICAL_CORES / 4);
    std::cout << "使用 " << optimal_threads << " 个工作线程" << std::endl;
    
    // 创建任务处理器
    TaskProcessor processor;
    std::vector<std::thread> threads;
    
    // 处理不同长度的序列
    const int max_possible_length = std::min({MAX_LENGTH, query_len, ref_len});
    
    for (int length = MIN_LENGTH; length <= max_possible_length; ++length) {
        if (length > query_len) continue;
        
        // 清理前一次迭代的数据
        processor.clear_positions();
        
        // 处理查询序列
        const int chunk_size = std::max(5000, query_len / (optimal_threads * 2));
        threads.clear();
        
        for (int i = 0; i <= query_len - length; i += chunk_size) {
            int end = std::min(i + chunk_size, query_len - length + 1);
            threads.emplace_back(&TaskProcessor::process_query_segment,
                               &processor, length, i, end);
        }
        
        for (auto& thread : threads) {
            thread.join();
        }
        
        // 处理参考序列
        const int ref_chunk_size = std::max(5000, ref_len / (optimal_threads * 2));
        threads.clear();
        
        for (int i = 0; i <= ref_len - length; i += ref_chunk_size) {
            int end = std::min(i + ref_chunk_size, ref_len - length + 1);
            threads.emplace_back(&TaskProcessor::process_reference_segment,
                               &processor, length, i, end);
        }
        
        for (auto& thread : threads) {
            thread.join();
        }
    }
    
    std::cout << std::endl << "所有任务处理完成，开始排序结果..." << std::endl;
    
    // 获取并处理结果
    auto& repeats = processor.get_results();
    std::sort(repeats.begin(), repeats.end());
    
    auto unique_end = std::unique(repeats.begin(), repeats.end(),
        [](const RepeatPattern& a, const RepeatPattern& b) {
            return a.position == b.position && 
                   a.is_reverse == b.is_reverse &&
                   a.length == b.length;
        });
    
    repeats.erase(unique_end, repeats.end());
    
    return repeats;
}

// 读取序列文件
std::string read_sequence(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("无法打开文件: " + filename);
    }
    
    std::string sequence;
    sequence.reserve(1024 * 1024); // 预分配1MB内存
    
    std::string line;
    while (std::getline(file, line)) {
        // 移除非DNA字符并转换为大写
        line.erase(
            std::remove_if(line.begin(), line.end(), 
                [](char c) { 
                    return !(c == 'A' || c == 'T' || c == 'G' || c == 'C' ||
                            c == 'a' || c == 't' || c == 'g' || c == 'c'); 
                }
            ),
            line.end()
        );
        
        std::transform(line.begin(), line.end(), line.begin(), ::toupper);
        sequence += line;
    }
    
    return sequence;
}

// 保存结果到文件
void save_repeats_to_file(const std::vector<RepeatPattern>& repeats, 
                         const std::string& output_file) {
    std::ofstream file(output_file);
    if (!file) {
        throw std::runtime_error("无法创建输出文件: " + output_file);
    }
    
    // 写入CSV头
    file << "参考位置,长度,重复次数,是否反向重复,原始序列,查询位置\n";
    
    // 写入重复片段信息
    for (const auto& repeat : repeats) {
        file << repeat.position << ","
             << repeat.length << ","
             << repeat.repeat_count << ","
             << (repeat.is_reverse ? "是" : "否") << ","
             << repeat.original_sequence << ","
             << repeat.query_position << "\n";
    }
    
    // 保存详细信息
    std::ofstream detail_file("repeat_details_stl.txt");
    if (!detail_file) {
        throw std::runtime_error("无法创建详细输出文件");
    }
    
    for (size_t i = 0; i < repeats.size(); ++i) {
        const auto& repeat = repeats[i];
        detail_file << "重复 #" << (i+1) << ":\n"
                   << "  参考位置: " << repeat.position << "\n"
                   << "  长度: " << repeat.length << "\n"
                   << "  重复次数: " << repeat.repeat_count << "\n"
                   << "  是否反向重复: " << (repeat.is_reverse ? "是" : "否") << "\n"
                   << "  原始序列: " << repeat.original_sequence << "\n"
                   << "  查询位置: " << repeat.query_position << "\n\n";
    }
}

int main(int argc, char* argv[]) {
    try {
        // 设置使用的线程数，选择合适的默认值
        // 对于数据密集应用，过多的线程可能导致资源竞争
        int num_threads = std::max(1, NUM_LOGICAL_CORES / 4);  
        std::cout << "请输入线程数 (默认" << num_threads << "): ";
        std::string input;
        std::getline(std::cin, input);
        
        if (!input.empty()) {
            try {
                num_threads = std::stoi(input);
                if (num_threads <= 0 || num_threads > NUM_LOGICAL_CORES) {
                    num_threads = std::max(1, NUM_LOGICAL_CORES / 4);
                }
            } catch (...) {
                // 保持默认值
            }
        }
        
        std::string query_file = "query.txt";
        std::string reference_file = "reference.txt";
        
        // 检查命令行参数
        if (argc >= 3) {
            reference_file = argv[1];
            query_file = argv[2];
        }
        
        std::cout << "读取查询序列: " << query_file << std::endl;
        std::string query = read_sequence(query_file);
        
        std::cout << "读取参考序列: " << reference_file << std::endl;
        std::string reference = read_sequence(reference_file);
        
        // 设置OpenMP线程数
        omp_set_num_threads(num_threads);
        
        // 开始计时 - 高精度计时器
        auto start = std::chrono::high_resolution_clock::now();
        
        // 查找重复
        auto repeats = find_repeats(query, reference);
        
        // 计算耗时（毫秒）
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        std::cout << "找到 " << repeats.size() << " 个重复片段，耗时: " 
                 << duration.count() << " 毫秒" << std::endl;
        
        // 输出重复片段详情
        for (const auto& repeat : repeats) {
            std::cout << "参考位置: " << repeat.position
                     << ", 长度: " << repeat.length
                     << ", 重复次数: " << repeat.repeat_count
                     << ", 是否反向重复: " << (repeat.is_reverse ? "是" : "否")
                     << ", 原始序列: " << repeat.original_sequence
                     << ", 查询位置: " << repeat.query_position << std::endl;
        }
        
        // 保存结果
        save_repeats_to_file(repeats, "repeat_results_stl.txt");
        
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}