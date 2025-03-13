#!/usr/bin/env python3
# DNA Repeat Finder
# 用于查找DNA序列中的重复片段
max_length = 120

def read_sequence(filename):
    """读取文件中的DNA序列"""
    with open(filename, 'r') as f:
        return f.read().strip().upper()

def get_reverse_complement(dna):
    """获取DNA序列的反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, 'N') for base in reversed(dna))

def build_lps(pattern):
    """构建 KMP 算法的最长相同前后缀数组"""
    length = len(pattern)
    lps = [0] * length
    len_prev = 0
    i = 1

    while i < length:
        if pattern[i] == pattern[len_prev]:
            len_prev += 1
            lps[i] = len_prev
            i += 1
        else:
            if len_prev != 0:
                len_prev = lps[len_prev - 1]
            else:
                lps[i] = 0
                i += 1
    return lps

def kmp_search(text, pattern):
    """使用 KMP 算法进行字符串匹配，返回所有匹配位置"""
    if not pattern or not text:
        return []

    matches = []
    n, m = len(text), len(pattern)
    lps = build_lps(pattern)
    
    i = j = 0
    while i < n:
        if pattern[j] == text[i]:
            i += 1
            j += 1
        if j == m:
            matches.append(i - j)
            j = lps[j - 1]
        elif i < n and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return matches

def count_occurrences(reference, segment):
    """计算片段在参考序列中的出现次数"""
    return len(kmp_search(reference, segment))

def find_repeats(query, reference, min_length=10):
    """查找重复片段，优化版本"""
    repeats = []
    ref_len = len(reference)
    query_len = len(query)
    
    print(f"搜索参数: min_length={min_length}, max_length={max_length}")
    
    # 使用滑动窗口和字典来跟踪片段位置
    window_positions = {}
    
    # 预处理query序列，记录所有可能片段的位置
    for length in range(min_length, min(max_length + 1, query_len - min_length + 1)):
        window_positions.clear()  # 清空字典以节省内存
        
        # 为每个长度构建一次窗口位置映射
        for i in range(query_len - length + 1):
            segment = query[i:i+length]
            if segment not in window_positions:
                window_positions[segment] = []
            window_positions[segment].append(i)
        
        # 检查reference中的每个片段
        for i in range(ref_len - length + 1):
            if i % 1000 == 0:
                print(f"处理进度: {i}/{ref_len} ({i/ref_len*100:.1f}%)")
            
            segment = reference[i:i+length]
            positions = window_positions.get(segment, [])
            
            # 检查正向重复
            if len(positions) >= 2:
                consecutive_positions = []
                current_group = [positions[0]]
                
                for j in range(1, len(positions)):
                    if positions[j] == current_group[-1] + length:
                        current_group.append(positions[j])
                    else:
                        if len(current_group) >= 2:
                            consecutive_positions.append(current_group[:])
                        current_group = [positions[j]]
                
                if len(current_group) >= 2:
                    consecutive_positions.append(current_group)
                
                for group in consecutive_positions:
                    repeats.append({
                        'position': i,
                        'length': length,
                        'repeat_count': len(group),
                        'is_reverse': False,
                        'original_sequence': segment,
                        'query_position': group[0]
                    })
            
            # 检查反向互补重复
            rev_comp = get_reverse_complement(segment)
            rev_positions = window_positions.get(rev_comp, [])
            
            if len(rev_positions) >= 2:
                consecutive_positions = []
                current_group = [rev_positions[0]]
                
                for j in range(1, len(rev_positions)):
                    if rev_positions[j] == current_group[-1] + length:
                        current_group.append(rev_positions[j])
                    else:
                        if len(current_group) >= 2:
                            consecutive_positions.append(current_group[:])
                        current_group = [rev_positions[j]]
                
                if len(current_group) >= 2:
                    consecutive_positions.append(current_group)
                
                for group in consecutive_positions:
                    repeats.append({
                        'position': i,
                        'length': length,
                        'repeat_count': len(group),
                        'is_reverse': True,
                        'original_sequence': segment,
                        'query_position': group[0]
                    })
    
    # 过滤嵌套重复并排序
    filtered_repeats = filter_nested_repeats(repeats)
    filtered_repeats.sort(key=lambda x: (x['length'] * x['repeat_count']), reverse=True)
    return filtered_repeats

def filter_nested_repeats(repeats):
    """
    过滤嵌套重复，保留每个位置最长的重复
    """
    if not repeats:
        return []
    
    # 按位置和反向标志分组
    position_groups = {}
    for repeat in repeats:
        key = (repeat['position'], repeat['is_reverse'])
        if key not in position_groups:
            position_groups[key] = []
        position_groups[key].append(repeat)
    
    # 从每组中选择最长的重复
    filtered_repeats = []
    for group in position_groups.values():
        best_repeat = max(group, key=lambda x: x['length'])
        filtered_repeats.append(best_repeat)
    
    return filtered_repeats

def save_repeats_to_file(repeats, output_file):
    """将重复片段信息保存到文件"""
    with open(output_file, 'w') as f:
        
        f.write("参考位置,长度,重复次数,是否反向重复,原始序列,查询位置\n")
        for repeat in repeats:
            repeat['position']+=repeat['length']
            is_reverse = "是" if repeat['is_reverse'] else "否"
            f.write(f"{repeat['position']},{repeat['length']},{repeat['repeat_count']},{is_reverse},{repeat['original_sequence'][:30]}...,{repeat['query_position']}\n")
    
    # 同时保存到详细文件
    with open('repeat_details.txt', 'w') as f:
        for i, repeat in enumerate(repeats):
            f.write(f"重复 #{i+1}:\n")
            f.write(f"  参考位置: {repeat['position']}\n")
            f.write(f"  长度: {repeat['length']}\n")
            f.write(f"  重复次数: {repeat['repeat_count']}\n")
            f.write(f"  是否反向重复: {'是' if repeat['is_reverse'] else '否'}\n")
            f.write(f"  原始序列: {repeat['original_sequence']}\n")
            f.write(f"  查询位置: {repeat['query_position']}\n\n")

def main():
    """主函数"""
    import time
    import sys
    
    # 处理命令行参数
    query_file = "query.txt"
    reference_file = "reference.txt"
    output_file = "repeat_results.txt"
    min_length = 10  # 默认最小重复片段长度
    
    # 检查是否有命令行参数
    if len(sys.argv) >= 3:
        reference_file = sys.argv[1]
        query_file = sys.argv[2]
    
    if len(sys.argv) >= 4:
        min_length = int(sys.argv[3])
    
    print(f"读取查询序列: {query_file}")
    query = read_sequence(query_file)
    print(f"查询序列长度: {len(query)}")
    
    print(f"读取参考序列: {reference_file}")
    reference = read_sequence(reference_file)
    print(f"参考序列长度: {len(reference)}")
    
    start_time = time.time()
    print(f"查找最小长度为 {min_length} 的重复片段...")
    repeats = find_repeats(query, reference, min_length)
    elapsed_time = time.time() - start_time
    
    print(f"找到 {len(repeats)} 个重复片段，耗时: {elapsed_time:.2f} 秒")
    
    print(f"保存结果到文件: {output_file}")
    save_repeats_to_file(repeats, output_file)
    
    # 输出前5个找到的重复片段
    print("\n前5个重复片段:")
    for i in range(min(5, len(repeats))):
        repeat = repeats[i]
        print(f"位置: {repeat['position']}, 长度: {repeat['length']}, 重复次数: {repeat['repeat_count']}, "
              f"是否反向: {'是' if repeat['is_reverse'] else '否'}, 序列: {repeat['original_sequence'][:20]}...")

if __name__ == "__main__":
    main()