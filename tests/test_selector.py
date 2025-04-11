import pytest
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.selector import FinalSelector
from sv_aligner.data_types import AlignmentSegment

class TestFinalSelector:
    def setup_method(self):
        """在每个测试方法前设置环境"""
        self.selector = FinalSelector()
    
    def test_select_non_overlapping_simple(self):
        """测试简单的非重叠段选择"""
        # 创建三个不重叠的段
        segments = [
            AlignmentSegment(q_name="test", q_len=100, q_st=0, q_en=10, r_name="ref", r_len=100, r_st=0, r_en=10, strand="+", score=20),
            AlignmentSegment(q_name="test", q_len=100, q_st=20, q_en=30, r_name="ref", r_len=100, r_st=20, r_en=30, strand="+", score=20),
            AlignmentSegment(q_name="test", q_len=100, q_st=40, q_en=50, r_name="ref", r_len=100, r_st=40, r_en=50, strand="+", score=20)
        ]
        
        # 选择段
        selected = self.selector.select_alignments(segments)
        
        # 验证结果 - 应该选择所有段
        assert len(selected) == 3
        assert all(seg in selected for seg in segments)
    
    def test_select_overlapping_segments(self):
        """测试处理重叠段"""
        # 创建三个段，其中第二个和第三个在查询序列上重叠
        segments = [
            AlignmentSegment(q_name="test", q_len=100, q_st=0, q_en=10, r_name="ref", r_len=100, r_st=0, r_en=10, strand="+", score=20),
            AlignmentSegment(q_name="test", q_len=100, q_st=20, q_en=40, r_name="ref", r_len=100, r_st=20, r_en=40, strand="+", score=30),
            AlignmentSegment(q_name="test", q_len=100, q_st=30, q_en=50, r_name="ref", r_len=100, r_st=30, r_en=50, strand="+", score=20)
        ]
        
        # 选择段
        selected = self.selector.select_alignments(segments)
        
        # 验证结果 - 应该选择第一个和第二个段（总分数50），而不是第一个和第三个段（总分数40）
        assert len(selected) == 2
        assert segments[0] in selected
        assert segments[1] in selected
        assert segments[2] not in selected
    
    def test_select_complex_overlapping(self):
        """测试复杂重叠情况"""
        # 创建一组复杂重叠的段
        segments = [
            # 第一个段（得分15）
            AlignmentSegment(q_name="test", q_len=100, q_st=0, q_en=20, r_name="ref", r_len=100, r_st=0, r_en=20, strand="+", score=15),
            # 第二个段（得分30）- 与第一个重叠
            AlignmentSegment(q_name="test", q_len=100, q_st=10, q_en=30, r_name="ref", r_len=100, r_st=30, r_en=50, strand="+", score=30),
            # 第三个段（得分20）- 与第四个重叠
            AlignmentSegment(q_name="test", q_len=100, q_st=40, q_en=60, r_name="ref", r_len=100, r_st=40, r_en=60, strand="+", score=20),
            # 第四个段（得分25）
            AlignmentSegment(q_name="test", q_len=100, q_st=50, q_en=70, r_name="ref", r_len=100, r_st=60, r_en=80, strand="+", score=25),
            # 第五个段（得分40）- 独立段
            AlignmentSegment(q_name="test", q_len=100, q_st=80, q_en=100, r_name="ref", r_len=100, r_st=80, r_en=100, strand="+", score=40)
        ]
        
        # 选择段
        selected = self.selector.select_alignments(segments)
        
        # 验证结果 - 最优解应该是第二个段(30) + 第四个段(25) + 第五个段(40)，总分数95
        assert len(selected) == 3
        assert segments[1] in selected
        assert segments[3] in selected
        assert segments[4] in selected
        
        # 验证总分数
        assert sum(segment.score for segment in selected) == 95
    
    def test_empty_segments(self):
        """测试空段列表"""
        selected = self.selector.select_alignments([])
        assert len(selected) == 0
    
    def test_single_segment(self):
        """测试单个段的选择"""
        segment = AlignmentSegment(q_name="test", q_len=100, q_st=0, q_en=10, r_name="ref", r_len=100, r_st=0, r_en=10, strand="+", score=20)
        selected = self.selector.select_alignments([segment])
        assert len(selected) == 1
        assert selected[0] == segment
