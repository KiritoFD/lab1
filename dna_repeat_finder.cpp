#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
#include <thread>
#include <future>
#include <locale>   // 添加 locale 头文件
#include <codecvt>  // 添加 codecvt 头文件
#include <chrono> // 添加计时头文件
#include <mutex>

// 重复信息的结构体
struct RepeatInfo {
    int position;
    int length;
    int count;
    bool is_reverse;
    std::string orig_seq;
    std::vector<std::string> repeat_examples;
};

// 生成DNA序列的反向互补序列
std::string get_reverse_complement(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.length());
    
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        switch (*it) {
            case 'A': result.push_back('T'); break;
            case 'T': result.push_back('A'); break;
            case 'G': result.push_back('C'); break;
            case 'C': result.push_back('G'); break;
            default: result.push_back('N'); // 处理非标准碱基
        }
    }
    
    return result;
}

// 构建相似度矩阵
std::vector<std::vector<int>> build_similarity_matrix(const std::string& reference, const std::string& query) {
    int n = reference.length();
    int m = query.length();
    std::vector<std::vector<int>> matrix(n, std::vector<int>(m, 0));
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            matrix[i][j] = (reference[i] == query[j]) ? 1 : -1;
        }
    }
    
    return matrix;
}

// 线程安全地添加重复信息
class ThreadSafeRepeatList {
public:
    void addRepeat(const RepeatInfo& repeat) {
        std::lock_guard<std::mutex> lock(mutex_);
        repeats_.push_back(repeat);
    }
    
    std::vector<RepeatInfo> getRepeats() const {
        std::lock_guard<std::mutex> lock(mutable_mutex_);
        return repeats_;
    }

private:
    std::vector<RepeatInfo> repeats_;
    mutable std::mutex mutable_mutex_;
    std::mutex mutex_;
};

// 并行查找重复
std::vector<RepeatInfo> find_repeats_parallel(const std::string& reference, const std::string& query, int num_threads) {
    ThreadSafeRepeatList safe_repeats;
    
    // 使用动态参数
    int min_length = std::max(5, static_cast<int>(reference.length() / 1000));
    int max_length = std::min(static_cast<int>(reference.length() / 10), 120);
    int step = std::max(1, min_length / 5);
    
    std::cout << "使用 " << num_threads << " 线程" << std::endl;
    
    // 将参考序列分成多个部分
    std::vector<std::pair<int, int>> ranges;
    int chunk_size = reference.length() / num_threads;
    for (int i = 0; i < num_threads; ++i) {
        int start = i * chunk_size;
        int end = (i == num_threads - 1) ? reference.length() : (i + 1) * chunk_size;
        ranges.push_back({start, end});
    }
    
    // 启动线程
    std::vector<std::future<void>> futures;
    for (const auto& range : ranges) {
        futures.push_back(std::async(std::launch::async, 
            [&](int start, int end) {
                for (int pos = start; pos <= end - min_length; ++pos) {
                    for (int length = min_length; length <= std::min(max_length, static_cast<int>(reference.length() - pos)); length += step) {
                        std::string segment = reference.substr(pos, length);
                        
                        // 检查正向重复
                        size_t start_idx = 0;
                        while (true) {
                            size_t next_idx = query.find(segment, start_idx);
                            if (next_idx == std::string::npos) {
                                break;
                            }
                            
                            // 检查后续是否有连续重复
                            size_t current_pos = next_idx + length;
                            int consecutive_count = 0;
                            
                            while (current_pos + length <= query.length()) {
                                if (query.substr(current_pos, length) == segment) {
                                    consecutive_count++;
                                    current_pos += length;
                                } else {
                                    break;
                                }
                            }
                            
                            if (consecutive_count > 0) {
                                // 检查是否是不同的序列
                                if (segment != query.substr(next_idx, length)) {
                                    RepeatInfo info;
                                    info.position = pos;
                                    info.length = length;
                                    info.count = consecutive_count;
                                    info.is_reverse = false;
                                    info.orig_seq = segment;
                                    
                                    // 添加重复实例
                                    for (int i = 0; i < consecutive_count; ++i) {
                                        info.repeat_examples.push_back(query.substr(next_idx + length * (i + 1), length));
                                    }
                                    
                                    safe_repeats.addRepeat(info);
                                }
                            }
                            
                            start_idx = next_idx + 1;
                            if (start_idx >= query.length()) {
                                break;
                            }
                        }
                        
                        // 检查反向互补重复
                        std::string rev_comp = get_reverse_complement(segment);
                        start_idx = 0;
                        
                        while (true) {
                            size_t next_idx = query.find(rev_comp, start_idx);
                            if (next_idx == std::string::npos) {
                                break;
                            }
                            
                            // 检查后续是否有连续重复
                            size_t current_pos = next_idx + length;
                            int consecutive_count = 0;
                            
                            while (current_pos + length <= query.length()) {
                                if (query.substr(current_pos, length) == rev_comp) {
                                    consecutive_count++;
                                    current_pos += length;
                                } else {
                                    break;
                                }
                            }
                            
                            // 反向互补必然是不同的
                            RepeatInfo info;
                            info.position = pos;
                            info.length = length;
                            info.count = std::max(1, consecutive_count);
                            info.is_reverse = true;
                            info.orig_seq = segment;
                            
                            // 添加重复实例
                            info.repeat_examples.push_back(query.substr(next_idx, length));
                            for (int i = 0; i < consecutive_count; ++i) {
                                info.repeat_examples.push_back(query.substr(next_idx + length * (i + 1), length));
                            }
                            
                            safe_repeats.addRepeat(info);
                            
                            start_idx = next_idx + 1;
                            if (start_idx >= query.length()) {
                                break;
                            }
                        }
                    }
                }
            }, range.first, range.second));
    }
    
    // 等待所有线程完成
    for (auto& future : futures) {
        future.get();
    }
    
    return safe_repeats.getRepeats();
}

// 寻找重复
std::vector<RepeatInfo> find_repeats(const std::string& reference, const std::string& query) {
    std::vector<RepeatInfo> repeats;
    
    // 使用动态参数
    int min_length = std::max(5, static_cast<int>(reference.length() / 1000));
    int max_length = std::min(static_cast<int>(reference.length() / 10), 120);
    int step = std::max(1, min_length / 5);
    
    // 对参考序列中的每个位置进行检查
    for (int pos = 0; pos <= reference.length() - min_length; ++pos) {
        // 对不同长度的片段进行检查
        for (int length = min_length; length <= std::min(max_length, static_cast<int>(reference.length() - pos)); length += step) {
            std::string segment = reference.substr(pos, length);
            
            // 检查正向重复
            size_t start_idx = 0;
            while (true) {
                size_t next_idx = query.find(segment, start_idx);
                if (next_idx == std::string::npos) {
                    break;
                }
                
                // 检查后续是否有连续重复
                size_t current_pos = next_idx + length;
                int consecutive_count = 0;
                
                while (current_pos + length <= query.length()) {
                    if (query.substr(current_pos, length) == segment) {
                        consecutive_count++;
                        current_pos += length;
                    } else {
                        break;
                    }
                }
                
                if (consecutive_count > 0) {
                    // 检查是否是不同的序列
                    if (segment != query.substr(next_idx, length)) {
                        RepeatInfo info;
                        info.position = pos;
                        info.length = length;
                        info.count = consecutive_count;
                        info.is_reverse = false;
                        info.orig_seq = segment;
                        
                        // 添加重复实例
                        for (int i = 0; i < consecutive_count; ++i) {
                            info.repeat_examples.push_back(query.substr(next_idx + length * (i + 1), length));
                        }
                        
                        repeats.push_back(info);
                    }
                }
                
                start_idx = next_idx + 1;
                if (start_idx >= query.length()) {
                    break;
                }
            }
            
            // 检查反向互补重复
            std::string rev_comp = get_reverse_complement(segment);
            start_idx = 0;
            
            while (true) {
                size_t next_idx = query.find(rev_comp, start_idx);
                if (next_idx == std::string::npos) {
                    break;
                }
                
                // 检查后续是否有连续重复
                size_t current_pos = next_idx + length;
                int consecutive_count = 0;
                
                while (current_pos + length <= query.length()) {
                    if (query.substr(current_pos, length) == rev_comp) {
                        consecutive_count++;
                        current_pos += length;
                    } else {
                        break;
                    }
                }
                
                // 反向互补必然是不同的
                RepeatInfo info;
                info.position = pos;
                info.length = length;
                info.count = std::max(1, consecutive_count);
                info.is_reverse = true;
                info.orig_seq = segment;
                
                // 添加重复实例
                info.repeat_examples.push_back(query.substr(next_idx, length));
                for (int i = 0; i < consecutive_count; ++i) {
                    info.repeat_examples.push_back(query.substr(next_idx + length * (i + 1), length));
                }
                
                repeats.push_back(info);
                
                start_idx = next_idx + 1;
                if (start_idx >= query.length()) {
                    break;
                }
            }
        }
    }
    
    return repeats;
}

// 过滤嵌套重复，只保留每个位置的最长重复
std::vector<RepeatInfo> filter_nested_repeats(const std::vector<RepeatInfo>& repeats) {
    // 按位置和方向分组
    std::map<std::pair<int, bool>, std::vector<RepeatInfo>> position_groups;
    
    for (const auto& repeat : repeats) {
        position_groups[{repeat.position, repeat.is_reverse}].push_back(repeat);
    }
    
    // 每组只保留最长的重复
    std::vector<RepeatInfo> filtered_repeats;
    for (auto& group : position_groups) {
        // 按长度降序排序
        std::sort(group.second.begin(), group.second.end(), 
              [](const RepeatInfo& a, const RepeatInfo& b) { return a.length > b.length; });
        
        // 保留最长的一个
        filtered_repeats.push_back(group.second[0]);
    }
    
    // 按位置排序
    std::sort(filtered_repeats.begin(), filtered_repeats.end(), 
          [](const RepeatInfo& a, const RepeatInfo& b) { return a.position < b.position; });
    
    return filtered_repeats;
}

// 使用动态规划寻找矩阵中的路径
std::vector<std::vector<std::tuple<int, int, int>>> find_paths_dp(
    const std::vector<std::vector<int>>& similarity_matrix,
    const std::string& reference,
    const std::string& query,
    int min_match_length = 10) {
    
    int n = similarity_matrix.size();
    int m = similarity_matrix[0].size();
    
    // DP表 - dp[i][j] = {score, length, prev_i, prev_j}
    std::vector<std::vector<std::tuple<int, int, int, int>>> dp(
        n, std::vector<std::tuple<int, int, int, int>>(m, {0, 0, -1, -1}));
    
    // 初始化第一个位置
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (similarity_matrix[i][j] == 1) {
                dp[i][j] = {1, 1, -1, -1};  // 没有前驱节点
            } else {
                dp[i][j] = {-1, 0, -1, -1};
            }
        }
    }
    
    // 填充DP表 - 对角线方向的匹配
    for (int i = 1; i < n; ++i) {
        for (int j = 1; j < m; ++j) {
            if (similarity_matrix[i][j] == 1) {
                auto [prev_score, prev_length, _, __] = dp[i-1][j-1];
                if (prev_score > 0) {  // 如果前一个位置有正得分，延续匹配
                    dp[i][j] = {prev_score + 1, prev_length + 1, i-1, j-1};
                }
            }
        }
    }
    
    // 找出所有有效的匹配路径（长度 >= min_match_length）
    std::vector<std::vector<std::tuple<int, int, int>>> all_paths;
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            auto [score, length, prev_i, prev_j] = dp[i][j];
            if (length >= min_match_length) {
                // 回溯构建路径
                std::vector<std::tuple<int, int, int>> path;
                int curr_i = i;
                int curr_j = j;
                
                while (curr_i >= 0 && curr_j >= 0) {
                    auto [curr_score, curr_length, prev_i, prev_j] = dp[curr_i][curr_j];
                    if (curr_score <= 0) break;
                    
                    path.emplace_back(curr_i, curr_j, 1);  // 1表示匹配方向
                    
                    if (prev_i < 0 || prev_j < 0) break;
                    curr_i = prev_i;
                    curr_j = prev_j;
                }
                
                // 路径是反向的，需要翻转
                std::reverse(path.begin(), path.end());
                
                // 只有长度达到要求时才保存
                if (path.size() >= min_match_length) {
                    // 检查是否有重复
                    int start_i = std::get<0>(path.front());
                    int start_j = std::get<1>(path.front());
                    int end_i = std::get<0>(path.back());
                    int end_j = std::get<1>(path.back());
                    int match_length = path.size();
                    
                    // 检查正向重复
                    std::vector<std::tuple<int, int, int>> repeat_paths;
                    int curr_j = end_j + 1;
                    
                    // 循环检查是否有连续重复
                    while (curr_j + match_length <= m) {
                        bool is_repeat = true;
                        std::vector<std::tuple<int, int, int>> temp_path;
                        
                        // 检查这段是否与原匹配段相同
                        for (int k = 0; k < match_length; ++k) {
                            if (curr_j + k >= m || query[start_j + k] != query[curr_j + k]) {
                                is_repeat = false;
                                break;
                            }
                            temp_path.emplace_back(start_i + k, curr_j + k, -1);  // -1表示重复
                        }
                        
                        if (is_repeat) {
                            repeat_paths.insert(repeat_paths.end(), temp_path.begin(), temp_path.end());
                            curr_j += match_length;
                        } else {
                            break;
                        }
                    }
                    
                    // 添加重复路径
                    if (!repeat_paths.empty()) {
                        std::vector<std::tuple<int, int, int>> combined_path = path;
                        combined_path.insert(combined_path.end(), repeat_paths.begin(), repeat_paths.end());
                        all_paths.push_back(combined_path);
                    }
                    
                    // 检查反向互补重复逻辑也可以类似添加
                }
            }
        }
    }
    
    return all_paths;
}

// 使用动态规划方法优化查找重复
std::vector<RepeatInfo> find_repeats_dp(const std::string& reference, const std::string& query) {
    // 构建相似度矩阵
    auto matrix = build_similarity_matrix(reference, query);
    
    // 使用DP查找路径
    auto paths = find_paths_dp(matrix, reference, query);
    
    // 从路径中提取重复信息
    std::vector<RepeatInfo> repeats;
    
    for (const auto& path : paths) {
        // 将路径分为匹配段和重复段
        std::vector<std::tuple<int, int>> match_segment;
        std::vector<std::tuple<int, int>> repeat_segment;
        
        for (const auto& [i, j, dir] : path) {
            if (dir == 1) {
                match_segment.push_back({i, j});
            } else if (dir == -1) {
                repeat_segment.push_back({i, j});
            }
        }
        
        if (!match_segment.empty()) {
            // 获取匹配段信息
            int start_i = std::get<0>(match_segment.front());
            int start_j = std::get<1>(match_segment.front());
            int end_i = std::get<0>(match_segment.back());
            int end_j = std::get<1>(match_segment.back());
            int segment_length = end_i - start_i + 1;
            
            // 计算重复次数
            int repeat_count = repeat_segment.size() / segment_length;
            
            if (repeat_count > 0) {
                RepeatInfo info;
                info.position = start_i;
                info.length = segment_length;
                info.count = repeat_count;
                info.is_reverse = false;
                info.orig_seq = reference.substr(start_i, segment_length);
                
                // 添加重复实例
                for (int i = 0; i < repeat_count; ++i) {
                    size_t idx = start_j + segment_length * (i + 1);
                    info.repeat_examples.push_back(query.substr(idx, segment_length));
                }
                
                repeats.push_back(info);
            }
        }
    }
    
    return repeats;
}

int main() {
    std::string reference;
    std::string query;
    
    // 尝试从文件读取序列
    try {
        std::ifstream ref_file("reference.txt");
        std::ifstream query_file("query.txt");
        
        if (ref_file && query_file) {
            std::getline(ref_file, reference);
            std::getline(query_file, query);
            std::cout << "成功从文件读取序列" << std::endl;
        } else {
            std::cout << "从文件读取失败，请输入序列" << std::endl;
            std::cout << "输入参考序列: ";
            std::cin >> reference;
            std::cout << "输入查询序列: ";
            std::cin >> query;
        }
    } catch (const std::exception& e) {
        std::cerr << "读取文件时出错: " << e.what() << std::endl;
        return 1;
    }
    
    // 确定线程数量 - 针对 AMD Ryzen 9 7940HX 进行优化
    // Ryzen 9 7940HX 有 8 核心 16 线程
    int num_threads = 16;
    
    // 记录开始时间
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // 查找重复 - 使用并行方法
    std::vector<RepeatInfo> repeats = find_repeats_parallel(reference, query, num_threads);
    
    // 记录结束时间
    auto end_time = std::chrono::high_resolution_clock::now();
    
    // 计算执行时间
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    // 使用动态规划方法查找重复
    std::vector<RepeatInfo> dp_repeats = find_repeats_dp(reference, query);
    
    // 合并两种方法的结果
    repeats.insert(repeats.end(), dp_repeats.begin(), dp_repeats.end());
    
    // 过滤嵌套重复
    std::vector<RepeatInfo> filtered_repeats = filter_nested_repeats(repeats);
    
    // 设置控制台输出编码为 UTF-8
    std::locale utf8_locale(std::locale(), new std::codecvt_utf8<wchar_t>);
    std::wcout.imbue(utf8_locale);
    
    // 输出结果
    std::wcout << L"\n找到的重复片段:" << std::endl;
    std::wcout << L"位置 | 长度 | 重复次数 | 是否反向重复 | 原始片段 | 重复实例" << std::endl;
    std::wcout << L"----------------------------------------------------------------------" << std::endl;
    
    for (const auto& repeat : filtered_repeats) {
        // 显示序列的前10个碱基，太长的用...表示
        std::wstring orig_display = std::wstring(repeat.orig_seq.begin(), repeat.orig_seq.end()).substr(0, 10);
        if (repeat.orig_seq.length() > 10) {
            orig_display += L"...";
        }
        
        // 显示重复实例
        std::wstring repeat_display;
        if (!repeat.repeat_examples.empty()) {
            repeat_display = std::wstring(repeat.repeat_examples[0].begin(), repeat.repeat_examples[0].end()).substr(0, 10);
            if (repeat.repeat_examples[0].length() > 10) {
                repeat_display += L"...";
            }
            
            if (repeat.repeat_examples.size() > 1) {
                repeat_display += L" (共" + std::to_wstring(repeat.repeat_examples.size()) + L"个实例)";
            }
        } else {
            repeat_display = L"未找到实例";
        }
        
        std::wprintf(L"%8d | %6d | %12d | %s | %s | %s\n", 
                    repeat.position, repeat.length, repeat.count, 
                    (repeat.is_reverse ? L"是" : L"否"), 
                    orig_display.c_str(), repeat_display.c_str());
    }
    
    std::wcout << L"\n共找到 " << filtered_repeats.size() << L" 个重复片段" << std::endl;
    
    
    
    // 保存结果到文件
    try {
        std::ofstream result_file("f:\\下载\\lab1\\repeat_results.txt");
        if (result_file.is_open()) { // 检查文件是否成功打开
            result_file << "找到的重复片段:\n";
            result_file << "位置 | 长度 | 重复次数 | 是否反向重复\n";
            result_file << "-------------------------------------------------\n";
            
            for (const auto& repeat : filtered_repeats) {
                result_file << std::setw(8) << repeat.position << " | "
                            << std::setw(6) << repeat.length << " | "
                            << std::setw(12) << repeat.count << " | "
                            << (repeat.is_reverse ? "是" : "否") << "\n";
            }
            
            result_file << "\n共找到 " << filtered_repeats.size() << " 个重复片段\n";
            std::cout << "\n基本重复片段信息已保存到 repeat_results.txt" << std::endl;
        } else {
            std::cerr << "无法打开文件 repeat_results.txt" << std::endl;
        }
        
        std::ofstream details_file("f:\\下载\\lab1\\repeat_details.txt");
        if (details_file) {
            details_file << "位置,长度,重复次数,是否反向重复,原始序列,重复实例\n";
            
            for (const auto& repeat : filtered_repeats) {
                details_file << repeat.position << ","
                            << repeat.length << ","
                            << repeat.count << ","
                            << (repeat.is_reverse ? "是" : "否") << ","
                            << repeat.orig_seq << ",";
                
                // 合并所有重复实例
                if (!repeat.repeat_examples.empty()) {
                    for (size_t i = 0; i < repeat.repeat_examples.size(); ++i) {
                        if (i > 0) details_file << ";";
                        details_file << repeat.repeat_examples[i];
                    }
                } else {
                    details_file << "未找到实例";
                }
                
                details_file << "\n";
            }
            
            std::cout << "详细的重复序列信息已保存到 repeat_details.txt" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "写入文件时出错: " << e.what() << std::endl;
    }
    // 输出执行时间
    std::wcout.imbue(std::locale("")); // 使用用户默认locale
    std::wcout << L"查找重复耗时: " << duration.count() << L" 毫秒" << std::endl;
    return 0;
}
