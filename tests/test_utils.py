import pytest
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.utils import reverse_complement, hash_kmer, get_minimizers

class TestUtils:
    """测试基本工具函数"""
    
    def test_reverse_complement(self):
        """测试序列反向互补功能"""
        assert reverse_complement("") == ""
        assert reverse_complement("A") == "T"
        assert reverse_complement("ACGT") == "ACGT"  # 自己的反向互补是自己
        assert reverse_complement("ACGTACGT") == "ACGTACGT"  # 回文序列
        assert reverse_complement("AATCG") == "CGATT"
        assert reverse_complement("acgt") == "ACGT"  # 测试小写输入
        assert reverse_complement("NNNN") == "NNNN"  # 测试N碱基
        assert reverse_complement("ACGTN") == "NACGT"  # 测试混合碱基
        assert reverse_complement("ATGC") == "GCAT"  # 基本测试

    def test_hash_kmer(self):
        """测试k-mer哈希功能"""
        # 相同k-mer应该产生相同的哈希值
        assert hash_kmer("ACGT") == hash_kmer("ACGT")
        # 不同k-mer应该产生不同的哈希值（在大多数情况下）
        assert hash_kmer("ACGT") != hash_kmer("ACGA")
        # 测试小写输入
        assert hash_kmer("acgt") == hash_kmer("ACGT")

    def test_get_minimizers(self):
        """测试最小化子提取功能"""
        # 测试基本的minimizer提取
        seq = "ACGTACGTACGT"
        k, w = 4, 3
        minimizers = get_minimizers(seq, k, w)
        assert len(minimizers) > 0  # 应该找到至少一个minimizer
        
        # 测试空序列
        assert get_minimizers("", k, w) == []
        
        # 测试序列长度小于k+w-1的情况
        assert get_minimizers("ACG", k, w) == []
        
        # 测试含有N的序列
        seq_with_n = "ACGTNACGT"
        minimizers_with_n = get_minimizers(seq_with_n, k, w)
        # 验证N被正确处理
        assert len(minimizers_with_n) > 0
        # N不应该出现在k-mer中
        for _, pos, _ in minimizers_with_n:
            k_mer = seq_with_n[pos:pos+k]
            assert 'N' not in k_mer
        
        # 测试不同w值的影响
        min_w1 = get_minimizers(seq, k, 1)
        min_w5 = get_minimizers(seq, k, 5)
        # 通常w越大，minimizers数量应该越少或相等
        assert len(min_w1) >= len(min_w5)
        
        # 测试具体minimizer的内容
        # 手动计算一个已知序列的minimizers进行比较
        simple_seq = "ACGTACGT"
        simple_k, simple_w = 3, 2
        simple_mins = get_minimizers(simple_seq, simple_k, simple_w)
        
        # 验证返回结果的格式
        for mini_hash, pos, is_forward in simple_mins:
            assert isinstance(mini_hash, int)
            assert isinstance(pos, int)
            assert isinstance(is_forward, bool)
            assert 0 <= pos <= len(simple_seq) - simple_k
            
        # 测试同一位置的k-mer在正反向链上的canonical选择
        palindrome = "ACGT"  # 自己是自己的反向互补
        pal_mins = get_minimizers(palindrome, 4, 1)
        # 应该只有一个minimizer
        assert len(pal_mins) == 1
        # 应该选择正向链(假设实现中默认选择正向链作为canonical k-mer)
        assert pal_mins[0][2] is True
