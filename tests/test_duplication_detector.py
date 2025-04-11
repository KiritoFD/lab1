import os
import sys
import pytest
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.data_types import AlignmentSegment
from sv_aligner.duplication_detector import (
    calculate_similarity,
    calculate_segment_similarity,
    detect_duplication,
    format_duplication_report,
    detect_and_report_duplications,
    DuplicationPattern
)

class TestDuplicationDetector:
    """测试拷贝数变异检测功能"""
    
    def setup_method(self):
        """每个测试前的设置"""
        self.alignments = self._create_test_alignments()
    
    def _create_test_alignments(self):
        """创建测试用的比对结果"""
        # 创建一个包含重复的例子：参考序列上的同一个区域在查询序列上出现两次
        ref_name = "reference"
        ref_len = 1000
        query_name = "query"
        query_len = 1200
        
        alignments = [
            # 第一个片段：查询序列中的第一个重复
            AlignmentSegment(
                q_name=query_name,
                q_len=query_len,
                q_st=100,
                q_en=200,  # 100bp长的片段
                r_name=ref_name,
                r_len=ref_len,
                r_st=500,
                r_en=600,  # 对应参考序列的500-600
                strand="+",
                score=100,
                edit_distance=0,
                cigar="100M"
            ),
            # 第二个片段：查询序列中的第二个重复，映射到同一参考位置
            AlignmentSegment(
                q_name=query_name,
                q_len=query_len,
                q_st=300,
                q_en=400,  # 另一个100bp的片段
                r_name=ref_name,
                r_len=ref_len,
                r_st=500,
                r_en=600,  # 映射到相同的参考区域
                strand="+",
                score=95,
                edit_distance=5,
                cigar="100M"
            ),
            # 第三个片段：无关的片段
            AlignmentSegment(
                q_name=query_name,
                q_len=query_len,
                q_st=500,
                q_en=700,
                r_name=ref_name,
                r_len=ref_len,
                r_st=700,
                r_en=900,
                strand="+",
                score=200,
                edit_distance=0,
                cigar="200M"
            ),
            # 第四个片段：第三个重复，但是反向互补的
            AlignmentSegment(
                q_name=query_name,
                q_len=query_len,
                q_st=800,
                q_en=900,
                r_name=ref_name,
                r_len=ref_len,
                r_st=500,
                r_en=600,
                strand="-",  # 反向链
                score=90,
                edit_distance=10,
                cigar="100M"
            )
        ]
        
        return alignments
    
    def test_calculate_similarity(self):
        """测试计算序列相似度功能"""
        assert calculate_similarity("ACGT", "ACGT") == 1.0  # 完全匹配
        assert calculate_similarity("ACGT", "ACGA") == 0.75  # 75%相似
        assert calculate_similarity("ACGT", "TGCA") == 0.0   # 完全不同
        assert calculate_similarity("", "") == 0.0        # 空字符串
        assert calculate_similarity("A", "AC") == 1.0     # 不同长度
    
    def test_calculate_segment_similarity(self):
        """测试计算比对片段相似度功能"""
        # 创建两个相似的片段
        seg1 = AlignmentSegment(
            q_name="query", q_len=100, q_st=0, q_en=50,
            r_name="ref", r_len=100, r_st=0, r_en=50,
            strand="+", score=100, edit_distance=0, cigar="50M"
        )
        
        seg2 = AlignmentSegment(
            q_name="query", q_len=100, q_st=50, q_en=100,
            r_name="ref", r_len=100, r_st=0, r_en=50,
            strand="+", score=90, edit_distance=5, cigar="50M"
        )
        
        # 计算相似度
        similarity = calculate_segment_similarity(seg1, seg2)
        assert similarity >= 0.0 and similarity <= 1.0
        
        # 编辑距离=0的片段与自己比较
        similarity_same = calculate_segment_similarity(seg1, seg1)
        assert similarity_same == 1.0
    
    def test_detect_duplication(self):
        """测试重复检测功能"""
        # 使用默认参数检测重复
        duplications = detect_duplication(self.alignments)
        
        # 应该至少找到一个重复
        assert len(duplications) >= 1
        
        # 验证找到的重复
        if len(duplications) > 0:
            dup = duplications[0]
            assert isinstance(dup, DuplicationPattern)
            
            # 检查属性
            assert dup.query_name == "query"
            # 验证它确实是从参考序列的同一区域映射而来
            assert dup.r_st1 == 500 and dup.r_en1 == 600
            assert dup.r_st2 == 500 and dup.r_en2 == 600
    
    def test_format_duplication_report(self):
        """测试格式化重复报告功能"""
        # 首先检测重复
        duplications = detect_duplication(self.alignments)
        
        if len(duplications) > 0:
            # 生成中文报告
            report_zh = format_duplication_report(duplications[0], use_chinese=True)
            assert "重复检测" in report_zh
            assert "查询序列" in report_zh
            
            # 生成英文报告
            report_en = format_duplication_report(duplications[0], use_chinese=False)
            assert "DUPLICATION DETECTED" in report_en
            assert "Query sequence" in report_en
    
    def test_detect_and_report_duplications(self):
        """测试检测并报告重复功能"""
        # 测试基本功能
        duplications = detect_and_report_duplications(
            self.alignments,
            min_length=50,
            min_similarity=0.8
        )
        
        assert len(duplications) >= 1
        
        # 测试更严格的相似度要求
        strict_dups = detect_and_report_duplications(
            self.alignments,
            min_length=50,
            min_similarity=0.99
        )
        
        # 应该找到更少的重复
        assert len(strict_dups) <= len(duplications)
