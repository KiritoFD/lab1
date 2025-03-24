#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

using namespace std;

// 重复片段的数据结构
struct RepeatPattern {
    int position;
    int length;
    int repeat_count;
    bool is_reverse;
    string original_sequence;
    int query_position;
};

class DNARepeatFinder {
public:
    // 超参数
    double similarity_threshold;
    int min_length;
    int max_length;
    int max_repeats;
    int prime;
    int mod;
    unordered_map<char, int> base_map;

    DNARepeatFinder(double similarity_threshold = 1.00,
                    int min_length = 50,
                    int max_length = 101,
                    int max_repeats = 1000,
                    int prime = 31,
                    int mod = (1 << 31) - 1,
                    unordered_map<char, int> base_map = {
                        {'A', 3}, {'C', 5}, {'G', 7}, {'T', 11}
                    }) :
            similarity_threshold(similarity_threshold),
            min_length(min_length),
            max_length(max_length),
            max_repeats(max_repeats),
            prime(prime),
            mod(mod),
            base_map(base_map) {
        _precompute_powers(max_length);
        // 预计算反向互补序列
        complement = {{'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}};
    }

private:
    vector<long long> powers;
    unordered_map<string, string> rev_comp_cache;
    unordered_map<char, char> complement;

    void _precompute_powers(int max_window) {
        powers.resize(max_window + 1);
        powers[0] = 1;
        for (int i = 1; i <= max_window; ++i) {
            powers[i] = (powers[i - 1] * prime) % mod;
        }
    }

    int rolling_hash(const string& sequence, int window_size) {
        long long hash_val = 0;
        for (int i = 0; i < window_size; ++i) {
            hash_val = (hash_val * prime + base_map[sequence[i]]) % mod;
        }
        return static_cast<int>(hash_val);
    }

    string get_reverse_complement(const string& sequence) {
        if (rev_comp_cache.count(sequence)) {
            return rev_comp_cache[sequence];
        }
        string rev_comp = "";
        for (int i = sequence.length() - 1; i >= 0; --i) {
            rev_comp += complement.count(sequence[i]) ? complement[sequence[i]] : 'N';
        }
        rev_comp_cache[sequence] = rev_comp;
        return rev_comp;
    }

    double calculate_similarity(const string& seq1, const string& seq2) {
        int matches = 0;
        int len = seq1.length();
        for (int i = 0; i < len; ++i) {
            if (seq1[i] == seq2[i]) {
                matches++;
            }
        }
        return static_cast<double>(matches) / len;
    }

    vector<int> find_consecutive_groups(const vector<int>& positions, int length) {
        if (positions.size() < 2) {
            return {};
        }

        vector<int> groups;
        int current_count = 1;

        for (size_t i = 1; i < positions.size(); ++i) {
            if (positions[i] == positions[i - 1] + length) {
                current_count++;
            } else {
                if (current_count >= 2) {
                    groups.push_back(current_count);
                }
                current_count = 1;
            }
        }

        if (current_count >= 2) {
            groups.push_back(current_count);
        }

        return groups;
    }

    vector<RepeatPattern> filter_nested_repeats(vector<RepeatPattern>& repeats) {
        if (repeats.size() <= 1) {
            return repeats;
        }

        vector<bool> to_remove(repeats.size(), false);
        vector<RepeatPattern> new_repeats;

        for (size_t i = 0; i < repeats.size(); ++i) {
            if (to_remove[i]) {
                continue;
            }

            for (size_t j = i + 1; j < repeats.size(); ++j) {
                if (to_remove[j]) {
                    continue;
                }

                if (repeats[i].position == repeats[j].position &&
                    repeats[i].is_reverse == repeats[j].is_reverse) {
                    // 保留较长的重复
                    if (repeats[i].length > repeats[j].length) {
                        to_remove[j] = true;
                    } else {
                        to_remove[i] = true;
                        break;
                    }
                }
            }
        }

        for (size_t i = 0; i < repeats.size(); ++i) {
            if (!to_remove[i]) {
                new_repeats.push_back(repeats[i]);
            }
        }

        return new_repeats;
    }

    unordered_map<int, vector<pair<string, int>>> build_sequence_hashmap(const string& sequence, int length) {
        unordered_map<int, vector<pair<string, int>>> window_positions;
        for (int i = 0; i <= (int)sequence.length() - length; ++i) {
            string segment = sequence.substr(i, length);
            int hash_val = rolling_hash(segment, length);
            window_positions[hash_val].push_back({segment, i});
        }
        return window_positions;
    }

public:
    vector<RepeatPattern> find_repeats(const string& query, const string& reference) {
        vector<RepeatPattern> repeats;
        int special_check_around = query.length() / 2;

        // 从环境变量获取特别关注区域位置
        const char* env_value = getenv("SPECIAL_CHECK_AREA");
        if (env_value) {
            try {
                special_check_around = stoi(env_value);
            } catch (const invalid_argument& e) {
                // 处理转换错误
            }
        }

        cout << "Special check area around position: " << special_check_around << endl;

        for (int length = min_length; length <= min(max_length + 1, (int)query.length()); ++length) {
            // 为当前长度构建查询序列的位置映射
            unordered_map<int, vector<pair<string, int>>> window_positions = build_sequence_hashmap(query, length);

            for (int i = 0; i <= (int)reference.length() - length; ++i) {
                // 对于非特殊区域限制检查长度范围
                if (abs(i - special_check_around) > 10 && length > min_length + 10) {
                    if (length > 100) {
                        break;
                    }
                }

                string segment = reference.substr(i, length);
                int segment_hash = rolling_hash(segment, length);

                // 检查正向重复
                vector<int> matches;
                if (similarity_threshold == 1.0) {
                    for (const auto& stored_seq_pos : window_positions[segment_hash]) {
                        if (stored_seq_pos.first == segment) {
                            matches.push_back(stored_seq_pos.second);
                        }
                    }
                } else {
                    for (const auto& stored_seq_pos : window_positions[segment_hash]) {
                        if (calculate_similarity(stored_seq_pos.first, segment) >= similarity_threshold) {
                            matches.push_back(stored_seq_pos.second);
                        }
                    }
                }

                if (matches.size() >= 2) {
                    vector<int> groups = find_consecutive_groups(matches, length);
                    for (int group_size : groups) {
                        if ((int)repeats.size() >= max_repeats) {
                            break;
                        }
                        repeats.push_back({i, length, group_size, false, segment, matches[0]});
                    }
                }

                // 检查反向互补重复
                string rev_comp = get_reverse_complement(segment);
                int rev_comp_hash = rolling_hash(rev_comp, length);

                vector<int> matches_rev;
                if (similarity_threshold == 1.0) {
                    for (const auto& stored_seq_pos : window_positions[rev_comp_hash]) {
                        if (stored_seq_pos.first == rev_comp) {
                            matches_rev.push_back(stored_seq_pos.second);
                        }
                    }
                } else {
                    for (const auto& stored_seq_pos : window_positions[rev_comp_hash]) {
                        if (calculate_similarity(stored_seq_pos.first, rev_comp) >= similarity_threshold) {
                            matches_rev.push_back(stored_seq_pos.second);
                        }
                    }
                }

                if (matches_rev.size() >= 2) {
                    vector<int> groups = find_consecutive_groups(matches_rev, length);
                    for (int group_size : groups) {
                        if ((int)repeats.size() >= max_repeats) {
                            break;
                        }
                        repeats.push_back({i, length, group_size, true, segment, matches_rev[0]});
                    }
                }
            }
        }

        // 过滤嵌套重复并排序
        repeats = filter_nested_repeats(repeats);
        sort(repeats.begin(), repeats.end(), [](const RepeatPattern& a, const RepeatPattern& b) {
            return a.length * a.repeat_count > b.length * b.repeat_count;
        });

        return repeats;
    }
};

string read_sequence(const string& filename) {
    int fd = open(filename.c_str(), O_RDONLY);
    if (fd == -1) {
        cerr << "Error opening file: " << filename << endl;
        return "";
    }

    long long file_size = lseek(fd, 0, SEEK_END);
    lseek(fd, 0, SEEK_SET);

    char* mapped_data = (char*)mmap(0, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped_data == MAP_FAILED) {
        cerr << "Error mapping file: " << filename << endl;
        close(fd);
        return "";
    }

    string sequence(mapped_data, file_size);

    if (munmap(mapped_data, file_size) == -1) {
        cerr << "Error unmapping file: " << filename << endl;
    }
    close(fd);

    string upper_sequence = "";
    for (char c : sequence) {
        if (isalpha(c)) {
            upper_sequence += toupper(c);
        }
    }

    return upper_sequence;
}

void save_repeats_to_file(const vector<RepeatPattern>& repeats, const string& base_filename) {
    string csv_filename = base_filename + "_results.csv";
    ofstream csv_file(csv_filename);
    csv_file << "Reference Position,Length,Repeat Count,Is Reverse Repeat,Original Sequence,Query Position\n";
    for (const auto& repeat : repeats) {
        csv_file << repeat.position + repeat.length << ","
                 << repeat.length << ","
                 << repeat.repeat_count << ","
                 << (repeat.is_reverse ? "Yes" : "No") << ","
                 << repeat.original_sequence << ","
                 << repeat.query_position << "\n";
    }
    csv_file.close();

    string details_filename = base_filename + "_details.txt";
    ofstream details_file(details_filename);
    for (size_t i = 0; i < repeats.size(); ++i) {
        details_file << "Repeat #" << i + 1 << ":\n";
        details_file << "  Reference Position: " << repeats[i].position << "\n";
        details_file << "  Length: " << repeats[i].length << "\n";
        details_file << "  Repeat Count: " << repeats[i].repeat_count << "\n";
        details_file << "  Is Reverse Repeat: " << (repeats[i].is_reverse ? "Yes" : "No") << "\n";
        details_file << "  Original Sequence: " << repeats[i].original_sequence << "\n";
        details_file << "  Query Position: " << repeats[i].query_position << "\n\n";
    }
    details_file.close();
}

int main(int argc, char* argv[]) {
    string query_file = "query.txt";
    string reference_file = "reference.txt";

    if (argc >= 3) {
        reference_file = argv[1];
        query_file = argv[2];
    }

    cout << "Reading query sequence: " << query_file << endl;
    string query = read_sequence(query_file);

    cout << "Reading reference sequence: " << reference_file << endl;
    string reference = read_sequence(reference_file);

    cout << "Query sequence length: " << query.length() << endl;
    cout << "Reference sequence length: " << reference.length() << endl;

    // 开始计时
    clock_t start_time = clock();

    // 创建DNA重复查找器实例并查找重复
    DNARepeatFinder finder;
    vector<RepeatPattern> repeats = finder.find_repeats(query, reference);

    // 计算耗时
    clock_t end_time = clock();
    double elapsed_time = double(end_time - start_time) / CLOCKS_PER_SEC * 1000;

    cout << "Found " << repeats.size() << " repeat fragments, elapsed time: " << elapsed_time << " ms" << endl;

    // 打印结果
    for (size_t i = 0; i < repeats.size(); ++i) {
        cout << "Repeat #" << i + 1 << ": Position " << repeats[i].position
             << ", Length " << repeats[i].length
             << ", Repeat Count " << repeats[i].repeat_count
             << ", Is Reverse Repeat " << (repeats[i].is_reverse ? "Yes" : "No") << endl;
    }

    // 保存结果
    save_repeats_to_file(repeats, "repeat_results_cpp");

    return 0;
}
