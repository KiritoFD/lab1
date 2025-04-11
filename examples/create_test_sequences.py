#!/usr/bin/env python3
"""
创建测试序列文件，包含各种结构变异
"""
import os
import random
import argparse

def generate_random_dna(length):
    """生成随机DNA序列"""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))

def reverse_complement(seq):
    """获取序列的反向互补"""
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(comp.get(base, 'N') for base in reversed(seq))

def save_sequence(seq, filename):
    """保存序列到文件"""
    with open(filename, 'w') as f:
        f.write(seq)
    print(f"Saved sequence to {filename} (length: {len(seq)})")

def create_reference_sequence(length=1000):
    """创建参考序列"""
    return generate_random_dna(length)

def create_query_with_sv(reference, options):
    """
    创建包含结构变异的查询序列
    
    Args:
        reference: 参考序列
        options: 控制创建什么类型的结构变异
    """
    ref_length = len(reference)
    query_parts = []
    
    # 起始匹配区域
    start_len = min(options.match_length, ref_length // 3)
    query_parts.append(reference[:start_len])
    
    current_pos = start_len
    
    # 添加缺失
    if options.deletion:
        deletion_size = min(options.deletion_size, ref_length // 5)
        current_pos += deletion_size  # 跳过参考序列中的这部分
    
    # 添加一段匹配
    middle_len = min(options.match_length, (ref_length - current_pos) // 2)
    if middle_len > 0:
        query_parts.append(reference[current_pos:current_pos+middle_len])
        current_pos += middle_len
    
    # 添加倒置
    if options.inversion:
        inversion_size = min(options.inversion_size, (ref_length - current_pos) // 4)
        if inversion_size > 0:
            inverted_seq = reverse_complement(reference[current_pos:current_pos+inversion_size])
            query_parts.append(inverted_seq)
            current_pos += inversion_size
    
    # 添加重复（拷贝数变异）
    if options.duplication:
        dup_start = random.randint(0, len(reference) // 2)
        dup_size = min(options.duplication_size, len(reference) // 4)
        dup_seq = reference[dup_start:dup_start+dup_size]
        query_parts.append(dup_seq)  # 添加重复区域
    
    # 添加易位
    if options.translocation:
        trans_start = max(0, ref_length - options.translocation_size - 100)
        trans_size = min(options.translocation_size, ref_length // 4)
        if trans_size > 0:
            trans_seq = reference[trans_start:trans_start+trans_size]
            query_parts.append(trans_seq)
    
    # 添加结尾匹配区域
    end_pos = min(current_pos + options.match_length, ref_length)
    if end_pos > current_pos:
        query_parts.append(reference[current_pos:end_pos])
    
    return ''.join(query_parts)

def main():
    parser = argparse.ArgumentParser(description='Create test sequences with structural variations')
    parser.add_argument('--ref-file', default='reference.txt', help='Output reference filename')
    parser.add_argument('--query-file', default='query.txt', help='Output query filename')
    parser.add_argument('--ref-length', type=int, default=1000, help='Reference sequence length')
    parser.add_argument('--match-length', type=int, default=100, help='Length of matching regions')
    parser.add_argument('--deletion', action='store_true', help='Include deletion')
    parser.add_argument('--deletion-size', type=int, default=100, help='Size of deletion')
    parser.add_argument('--inversion', action='store_true', help='Include inversion')
    parser.add_argument('--inversion-size', type=int, default=100, help='Size of inversion')
    parser.add_argument('--duplication', action='store_true', help='Include duplication')
    parser.add_argument('--duplication-size', type=int, default=100, help='Size of duplication')
    parser.add_argument('--translocation', action='store_true', help='Include translocation')
    parser.add_argument('--translocation-size', type=int, default=100, help='Size of translocation')
    
    args = parser.parse_args()
    
    # 创建参考序列
    reference = create_reference_sequence(args.ref_length)
    save_sequence(reference, args.ref_file)
    
    # 创建包含结构变异的查询序列
    query = create_query_with_sv(reference, args)
    save_sequence(query, args.query_file)

if __name__ == "__main__":
    main()
