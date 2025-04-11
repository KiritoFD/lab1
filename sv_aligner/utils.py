import sys
from typing import List, Tuple, Dict, Iterator

def reverse_complement(seq: str) -> str:
    """计算DNA序列的反向互补"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    # 处理可能存在的小写字母
    seq_upper = seq.upper()
    return "".join(complement.get(base, base) for base in reversed(seq_upper))

def hash_kmer(kmer: str) -> int:
    """计算k-mer的哈希值 (使用简单的多项式哈希)"""
    val = 0
    base = 4
    for char in kmer.upper():
        # 使用简单的base-4表示法
        idx = 'ACGT'.find(char)
        if idx != -1:  # 只处理ACGT碱基
            val = val * base + idx
        else:
            val = val * base  # 对于非ACGT碱基，视为A (值为0)
    return val

def get_minimizers(seq: str, k: int, w: int) -> List[Tuple[int, int, bool]]:
    """
    从序列中提取minimizers。
    
    Args:
        seq: DNA序列
        k: k-mer大小
        w: 窗口大小
    
    Returns:
        List of tuples: (minimizer_hash, position, is_forward_strand)
        'is_forward_strand'表示产生canonical minimizer的k-mer是否在正向链。
    """
    if len(seq) < k + w - 1:
        return []

    minimizers = []
    
    # 滑动窗口处理序列
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k].upper()
        # 跳过含有N的k-mer
        if 'N' in kmer:
            continue

        kmer_rc = reverse_complement(kmer)

        # 使用canonical k-mer (正向和反向互补中字典序较小的一个)
        if kmer <= kmer_rc:
            canonical_kmer = kmer
            is_forward = True
        else:
            canonical_kmer = kmer_rc
            is_forward = False

        current_hash = hash_kmer(canonical_kmer)

        # 标准minimizer逻辑
        # 在窗口[i-w+1, i]范围内检查当前位置是否是最小值
        window_start = max(0, i - w + 1)
        min_hash_in_window = float('inf')
        pos_of_min_hash = -1
        is_fwd_of_min = True

        for j in range(window_start, i + 1):
            win_kmer = seq[j:j+k].upper()
            if 'N' in win_kmer:
                continue
            win_kmer_rc = reverse_complement(win_kmer)
            win_canonical = win_kmer if win_kmer <= win_kmer_rc else win_kmer_rc
            win_hash = hash_kmer(win_canonical)
            win_is_fwd = win_kmer <= win_kmer_rc

            if win_hash < min_hash_in_window:
                min_hash_in_window = win_hash
                pos_of_min_hash = j
                is_fwd_of_min = win_is_fwd
            # 如果哈希值相等，保留左边的（即较早出现的）

        # 如果此窗口的minimizer正好是当前位置，且之前没有相同哈希/位置添加过
        if pos_of_min_hash == i:
            if not minimizers or minimizers[-1][:2] != (min_hash_in_window, pos_of_min_hash):
                minimizers.append((min_hash_in_window, pos_of_min_hash, is_fwd_of_min))

    return minimizers
