import pytest
from unittest.mock import MagicMock, patch
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.refiner import Refiner
from sv_aligner.data_types import Seed, Chain, AlignmentSegment

class TestRefiner:
    def setup_method(self):
        """在每个测试方法前设置环境"""
        # 创建参数字典
        self.params = {
            'match_score': 2,
            'mismatch_penalty': -3,
            'gap_open_penalty': -5,
            'gap_extend_penalty': -2,
            'k': 4,
            'min_segment_length': 5,
            'min_alignment_score': 10
        }
        
        # 创建mock对象
        self.mock_aligner = MagicMock()
        
        # 准备参考序列
        self.ref_sequences = {
            0: "ACGTACGTACGTACGTACGT",  # ref_id 0
            1: "TGCATGCATGCATGCATGCA"   # ref_id 1
        }
        
        # 准备参考序列信息
        self.ref_info = {
            0: ("ref1", 20),  # ref1长度为20
            1: ("ref2", 20)   # ref2长度为20
        }
        
        # 创建refiner实例
        self.refiner = Refiner(self.mock_aligner, self.ref_sequences, self.ref_info, self.params)
    
    def test_refine_simple_chain(self):
        """测试简单链的精细化"""
        # 配置mock aligner返回值
        self.mock_aligner.align.return_value = (16, "14M", 0, 14, 0, 14)
        
        # 创建一个简单的正向链
        seeds = [
            Seed(q_pos=0, ref_id=0, r_pos=0, strand=True),
            Seed(q_pos=5, ref_id=0, r_pos=5, strand=True),
            Seed(q_pos=10, ref_id=0, r_pos=10, strand=True)
        ]
        
        chain = Chain(
            seeds=seeds,
            q_start=0,
            q_end=14,
            r_start=0,
            r_end=14,
            strand=True,
            score=15.0
        )
        
        # 创建一个查询序列，与参考序列的对应部分完全匹配
        query = self.ref_sequences[0][:14]
        
        # 精细化链
        segments = self.refiner.refine_chains([chain], query, "test_query")
        
        # 验证aligner.align被调用
        self.mock_aligner.align.assert_called_once()
        
        # 验证结果
        assert len(segments) == 1
        segment = segments[0]
        
        # 验证段的基本属性
        assert segment.q_name == "test_query"
        assert segment.q_st == 0
        assert segment.q_en == 14
        assert segment.r_name == "ref1"
        assert segment.r_st == 0
        assert segment.r_en == 14
        assert segment.strand == '+'
        assert segment.score == 16
    
    def test_refine_reverse_chain(self):
        """测试反向链的精细化"""
        # 配置mock aligner返回值
        self.mock_aligner.align.return_value = (16, "14M", 0, 14, 0, 14)
        
        # 创建一个反向链
        seeds = [
            Seed(q_pos=0, ref_id=0, r_pos=10, strand=False),
            Seed(q_pos=5, ref_id=0, r_pos=5, strand=False),
            Seed(q_pos=10, ref_id=0, r_pos=0, strand=False)
        ]
        
        chain = Chain(
            seeds=seeds,
            q_start=0,
            q_end=14,
            r_start=0,
            r_end=14,
            strand=False,
            score=15.0
        )
        
        # 创建一个查询序列
        query = "ACGTACGTACGTAC"  # 任意序列，因为我们mock了aligner
        
        # 精细化链
        segments = self.refiner.refine_chains([chain], query, "test_query")
        
        # 验证aligner.align被调用
        self.mock_aligner.align.assert_called_once()
        
        # 验证结果
        assert len(segments) == 1
        segment = segments[0]
        
        # 验证段的基本属性
        assert segment.q_name == "test_query"
        assert segment.strand == '-'
        assert segment.score == 16
    
    def test_refine_filter_by_length_and_score(self):
        """测试根据长度和分数过滤段"""
        # 配置mock aligner返回值
        self.mock_aligner.align.return_value = (5, "4M", 0, 4, 0, 4)  # 低分数、短长度
        
        # 创建一个应该被过滤掉的短链
        seeds = [
            Seed(q_pos=0, ref_id=0, r_pos=0, strand=True),
            Seed(q_pos=2, ref_id=0, r_pos=2, strand=True)
        ]
        
        short_chain = Chain(
            seeds=seeds,
            q_start=0,
            q_end=4,
            r_start=0,
            r_end=4,
            strand=True,
            score=5.0
        )
        
        # 创建查询序列
        query = "ACGT"
        
        # 设置较高的过滤阈值
        strict_params = dict(self.params)
        strict_params['min_segment_length'] = 10
        strict_params['min_alignment_score'] = 20
        
        strict_refiner = Refiner(self.mock_aligner, self.ref_sequences, self.ref_info, strict_params)
        
        # 精细化链
        segments = strict_refiner.refine_chains([short_chain], query, "test_query")
        
        # 验证结果 - 应该没有通过过滤
        assert len(segments) == 0
    
    def test_calculate_edit_distance_from_cigar(self):
        """测试从CIGAR字符串计算编辑距离"""
        # 测试各种CIGAR字符串
        assert self.refiner._calculate_edit_distance_from_cigar("10M") == 0
        assert self.refiner._calculate_edit_distance_from_cigar("5M2I3M") == 2
        assert self.refiner._calculate_edit_distance_from_cigar("5M2D3M") == 2
        assert self.refiner._calculate_edit_distance_from_cigar("5M2I3M1D2M") == 3
        assert self.refiner._calculate_edit_distance_from_cigar("") == 0
        assert self.refiner._calculate_edit_distance_from_cigar("5=2X3=") == 0  # 只计算I和D
