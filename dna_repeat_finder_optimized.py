#!/usr/bin/env python3
# DNA Repeat Finder - 优化版本
# 时间复杂度: O(|R| × |Q|) 代替原来的 O(|R|² × |Q|)
import time
import sys
from collections import defaultdict

max_length = 120  # 最大重复片段长度
min_length = 10   # 默认最小重复片段长度

def read_sequence(filename):
    """读取文件中的DNA序列"""
    with open(filename, 'r') as f:
        return f.read().strip().upper()  # 转换为大写以统一处理

def get_reverse_complement(dna):
    """获取DNA序列的反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, 'N') for base in reversed(dna))

def create_reference_index(reference, min_length):
    """
    为参考序列创建索引，记录每个片段的出现位置和次数
    这个预处理步骤让我们可以在O(1)时间内检查参考序列中的唯一性
    
    返回:
        segment_positions: 每个片段的位置列表
        segment_counts: 每个片段的出现次数
    """
    segment_positions = defaultdict(list)
    segment_counts = {}
    ref_len = len(reference)
    
    print("预处理参考序列...")
    # 为每个可能的min_length长度的片段建立索引
    for pos in range(ref_len - min_length + 1):
        if pos % 10000 == 0:
            print(f"索引处理进度: {pos}/{ref_len} ({pos/ref_len*100:.1f}%)")
        
        segment = reference[pos:pos+min_length]
        segment_positions[segment].append(pos)
    
    # 计算每个片段的出现次数
    for segment, positions in segment_positions.items():
        segment_counts[segment] = len(positions)
    
    return segment_positions, segment_counts

def find_connected_seeds(segment_positions, segment_counts, query, reference, min_length):
    """
    使用种子扩展法寻找潜在的重复片段
    时间复杂度: O(|R| + |Q|)
    
    参数:
        segment_positions: 每个片段的位置字典
        segment_counts: 每个片段的出现次数字典
        query: 查询DNA序列
        reference: 参考DNA序列
        min_length: 最小重复片段长度
        
    返回:
        potential_repeats: 潜在重复片段的列表，包含位置和长度信息
    """
    potential_repeats = []
    query_len = len(query)
    
    # 在查询序列中滑动窗口
    for i in range(query_len - min_length + 1):
        if i % 10000 == 0:
            print(f"查询处理进度: {i}/{query_len} ({i/query_len*100:.1f}%)")
        
        # 获取查询序列中的当前片段
        query_segment = query[i:i+min_length]
        
        # 检查这个片段是否在参考序列中且只出现一次
        if query_segment in segment_counts and segment_counts[query_segment] == 1:
            # 找到参考序列中的位置
            ref_pos = segment_positions[query_segment][0]
            
            # 尝试向右扩展以找到最长匹配
            extend_length = extend_match(query, reference, i, ref_pos, min_length)
            if extend_length >= min_length:
                potential_repeats.append({
                    'position': ref_pos,
                    'length': extend_length,
                    'query_position': i
                })
                
        # 同样检查反向互补
        rev_comp = get_reverse_complement(query_segment)
        if rev_comp in segment_counts and segment_counts[rev_comp] == 1:
            ref_pos = segment_positions[rev_comp][0]
            # 对于反向互补，扩展方式稍有不同
            extend_length = extend_match_reverse_complement(
                query, reference, i, ref_pos, min_length
            )
            if extend_length >= min_length:
                potential_repeats.append({
                    'position': ref_pos,
                    'length': extend_length,
                    'query_position': i,
                    'is_reverse': True
                })
    
    return potential_repeats

def extend_match(query, reference, query_pos, ref_pos, seed_length):
    """
    从种子匹配开始向两边扩展，找到最长匹配
    """
    query_len = len(query)
    ref_len = len(reference)
    
    # 初始长度为种子长度
    length = seed_length
    
    # 向右扩展
    q_right = query_pos + seed_length
    r_right = ref_pos + seed_length
    
    while q_right < query_len and r_right < ref_len and query[q_right] == reference[r_right]:
        length += 1
        q_right += 1
        r_right += 1
    
    # 向左扩展
    q_left = query_pos - 1
    r_left = ref_pos - 1
    
    while q_left >= 0 and r_left >= 0 and query[q_left] == reference[r_left]:
        length += 1
        q_left -= 1
        r_left -= 1
    
    return length

def extend_match_reverse_complement(query, reference, query_pos, ref_pos, seed_length):
    """
    为反向互补匹配特别设计的扩展函数
    """
    query_len = len(query)
    ref_len = len(reference)
    
    # 反向互补的种子长度
    length = seed_length
    
    # 向右扩展查询序列，向左扩展参考序列（因为是反向互补）
    q_right = query_pos + seed_length
    r_left = ref_pos - 1
    
    while q_right < query_len and r_left >= 0:
        # 获取参考序列中的碱基的互补
        ref_base = reference[r_left]
        comp_base = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}.get(ref_base, 'N')
        
        if query[q_right] == comp_base:
            length += 1
            q_right += 1
            r_left -= 1
        else:
            break
    
    # 向左扩展查询序列，向右扩展参考序列
    q_left = query_pos - 1
    r_right = ref_pos + seed_length
    
    while q_left >= 0 and r_right < ref_len:
        ref_base = reference[r_right]
        comp_base = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}.get(ref_base, 'N')
        
        if query[q_left] == comp_base:
            length += 1
            q_left -= 1
            r_right += 1
        else:
            break
    
    return length

def check_consecutive_repeats(potential_repeats, query, reference):
    """
    从潜在重复片段中检查连续重复
    时间复杂度: O(|Q|)
    """
    repeats = []
    query_len = len(query)
    
    for repeat in potential_repeats:
        position = repeat['position']
        length = repeat['length']
        query_position = repeat['query_position']
        is_reverse = repeat.get('is_reverse', False)
        
        if is_reverse:
            segment = get_reverse_complement(reference[position:position+length])
        else:
            segment = reference[position:position+length]
        
        # 检查是否有连续重复
        count = 1
        current_pos = query_position + length
        
        while current_pos + length <= query_len:
            if query[current_pos:current_pos+length] == segment:
                count += 1
                current_pos += length
            else:
                break
        
        if count > 1:  # 确保至少重复两次
            repeats.append({
                'position': position,
                'length': length,
                'repeat_count': count,
                'is_reverse': is_reverse,
                'original_sequence': reference[position:position+length],
                'query_position': query_position
            })
    
    return repeats

def filter_nested_repeats(repeats):
    """
    过滤嵌套重复，保留每个位置最长的重复
    """
    if not repeats:
        return []
    
    # 按位置和反向标志分组
    position_groups = {}
    for repeat in repeats:
        key = (repeat['position'], repeat.get('is_reverse', False))
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
            # 调整位置以保持与原始算法兼容
            position = repeat['position'] + repeat['length']
            is_reverse = "是" if repeat.get('is_reverse', False) else "否"
            f.write(f"{position},{repeat['length']},{repeat['repeat_count']},{is_reverse},{repeat['original_sequence'][:30]}...,{repeat['query_position']}\n")
    
    # 同时保存到详细文件
    with open('repeat_details.txt', 'w') as f:
        for i, repeat in enumerate(repeats):
            f.write(f"重复 #{i+1}:\n")
            f.write(f"  参考位置: {repeat['position']}\n")
            f.write(f"  长度: {repeat['length']}\n")
            f.write(f"  重复次数: {repeat['repeat_count']}\n")
            f.write(f"  是否反向重复: {'是' if repeat.get('is_reverse', False) else '否'}\n")
            f.write(f"  原始序列: {repeat['original_sequence']}\n")
            f.write(f"  查询位置: {repeat['query_position']}\n\n")

def find_repeats_optimized(query, reference, min_length=10):
    """
    优化版DNA重复查找算法，时间复杂度O(|R| × |Q|)
    
    参数:
        query: 查询DNA序列
        reference: 参考DNA序列
        min_length: 最小重复片段长度
        
    返回:
        重复片段列表，每个重复包含位置、长度、重复次数等信息
    """
    print(f"优化算法启动: min_length={min_length}, max_length={max_length}")
    
    # 步骤1: 为参考序列创建索引 - O(|R|)
    segment_positions, segment_counts = create_reference_index(reference, min_length)
    
    # 步骤2: 使用种子扩展查找潜在重复 - O(|Q|)
    potential_repeats = find_connected_seeds(
        segment_positions, segment_counts, query, reference, min_length
    )
    
    print(f"找到 {len(potential_repeats)} 个潜在重复片段")
    
    # 步骤3: 检查连续重复 - O(|Q|)
    repeats = check_consecutive_repeats(potential_repeats, query, reference)
    
    print(f"找到 {len(repeats)} 个确认重复片段")
    
    # 步骤4: 过滤和排序 - O(N log N)，N是重复片段数
    filtered_repeats = filter_nested_repeats(repeats)
    filtered_repeats.sort(
        key=lambda x: (x['length'] * x['repeat_count']), 
        reverse=True
    )
    
    return filtered_repeats

def main():
    """主函数"""
    # 处理命令行参数
    query_file = "query.txt"
    reference_file = "reference.txt"
    output_file = "repeat_results_optimized.txt"
    global min_length
    
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
    print(f"使用优化算法查找最小长度为 {min_length} 的重复片段...")
    repeats = find_repeats_optimized(query, reference, min_length)
    elapsed_time = time.time() - start_time
    
    print(f"找到 {len(repeats)} 个重复片段，耗时: {elapsed_time:.2f} 秒")
    
    print(f"保存结果到文件: {output_file}")
    save_repeats_to_file(repeats, output_file)
    
    # 输出前5个找到的重复片段
    print("\n前5个重复片段:")
    for i in range(min(5, len(repeats))):
        repeat = repeats[i]
        print(f"位置: {repeat['position']}, 长度: {repeat['length']}, 重复次数: {repeat['repeat_count']}, "
              f"是否反向: {'是' if repeat.get('is_reverse', False) else '否'}, 序列: {repeat['original_sequence'][:20]}...")

if __name__ == "__main__":
    main()