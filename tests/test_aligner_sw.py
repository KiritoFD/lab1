import pytest
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.aligner_sw import AlignerSW

class TestAlignerSW:
    def setup_method(self):
        """在每个测试方法前设置环境"""
        # 创建参数字典
        self.params = {
            'match_score': 2,
            'mismatch_penalty': -3,
            'gap_open_penalty': -5,
            'gap_extend_penalty': -2
        }
        
        # 创建aligner实例
        self.aligner = AlignerSW(self.params)
    
    def test_align_exact_match(self):
        """测试完全匹配的序列对齐"""
        # 两个完全匹配的序列
        query = "ACGTACGT"
        ref = "ACGTACGT"
        
        result = self.aligner.align(query, ref)
        
        # 判断返回值是否是有效的元组
        assert result is not None
        score, cigar, q_st, q_en, r_st, r_en = result
        
        # 校验结果
        assert q_st == 0
        assert q_en == 8
        assert r_st == 0
        assert r_en == 8
        assert score > 0
        # 完全匹配应该是8M或8=
        assert "8M" in cigar or "8=" in cigar
    
    def test_align_with_mismatches(self):
        """测试含有错配的序列对齐"""
        # 有2个错配的序列
        query = "ACGTAGGT"  # 第5个和第6个位置有错配
        ref = "ACGTACGT"
        
        result = self.aligner.align(query, ref)
        
        # 判断返回值是否是有效的元组
        assert result is not None
        score, cigar, q_st, q_en, r_st, r_en = result
        
        # 校验结果
        assert q_en - q_st == 8  # 应该覆盖整个查询序列
        assert r_en - r_st == 8  # 应该覆盖整个参考序列
        
        # 错配应该有CIGAR中的X操作或者在M操作中
        # 但具体CIGAR字符串格式取决于Smith-Waterman实现
    
    def test_align_with_insertion(self):
        """测试含有插入的序列对齐"""
        # 查询序列中有一个2bp插入
        query = "ACGTACGTAC"  # 末尾多了AC
        ref = "ACGTACGT"
        
        result = self.aligner.align(query, ref)
        
        # 判断返回值是否是有效的元组
        assert result is not None
        score, cigar, q_st, q_en, r_st, r_en = result
        
        # 校验结果
        assert q_en - q_st == 10  # 完整查询序列长度
        assert r_en - r_st == 8   # 完整参考序列长度
        
        # CIGAR应该包含插入或有M操作
        assert "I" in cigar or "M" in cigar
    
    def test_align_with_deletion(self):
        """测试含有缺失的序列对齐"""
        # 查询序列中有一个2bp缺失
        query = "ACGTGT"
        ref = "ACGTACGT"
        
        result = self.aligner.align(query, ref)
        
        # 判断返回值是否是有效的元组
        assert result is not None
        score, cigar, q_st, q_en, r_st, r_en = result
        
        # 校验结果
        assert q_en - q_st == 6   # 完整查询序列长度
        
        # CIGAR应该包含缺失或有M操作
        assert "D" in cigar or "M" in cigar
    
    def test_align_empty_sequences(self):
        """测试空序列的处理"""
        # 两个空序列
        result = self.aligner.align("", "")
        assert result is None
        
        # 一个空序列
        result = self.aligner.align("ACGT", "")
        assert result is None
        result = self.aligner.align("", "ACGT")
        assert result is None
