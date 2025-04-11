import os
import pytest
import tempfile
from unittest.mock import patch
from sv_aligner.orchestrator import Orchestrator
from sv_aligner.data_types import AlignmentSegment

class TestOrchestrator:
    def setup_method(self):
        """在每个测试方法前设置环境"""
        # 创建临时目录
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # 创建参数字典
        self.params = {
            'k': 4,
            'w': 2,
            'match_score': 2,
            'mismatch_penalty': -3,
            'gap_open_penalty': -5,
            'gap_extend_penalty': -2,
            'min_segment_length': 5,
            'min_alignment_score': 10,
            'min_chain_score': 10,
            'gap_penalty_factor': 0.01,
            'dist_diff_penalty_factor': 0.05,
            'max_gap': 1000,
            'max_dist_diff': 500,
            'min_query_len': 10
        }
        
        # 创建测试文件
        self.ref_fasta = os.path.join(self.temp_dir.name, "ref.fa")
        self.query_fasta = os.path.join(self.temp_dir.name, "query.fa")
        
        # 写入测试参考序列
        with open(self.ref_fasta, "w") as f:
            f.write(">ref1\n")
            f.write("ACGTACGTACGTACGTACGT\n")
            f.write(">ref2\n")
            f.write("TGCATGCATGCATGCATGCA\n")
        
        # 写入测试查询序列
        with open(self.query_fasta, "w") as f:
            f.write(">query1\n")
            f.write("ACGTACGTACGT\n")
            f.write(">query2\n")
            f.write("TGCATGCATGCA\n")
    
    def teardown_method(self):
        """在每个测试方法后清理环境"""
        self.temp_dir.cleanup()
    
    def test_load_reference_and_index(self):
        """测试参考序列加载和索引构建"""
        orchestrator = Orchestrator(self.params)
        
        # 加载参考序列和构建索引
        orchestrator.load_reference_and_index(self.ref_fasta)
        
        # 验证参考序列和索引已加载
        assert orchestrator.ref_index is not None
        assert orchestrator.ref_info is not None
        assert orchestrator.ref_sequences is not None
        
        # 验证参考序列数量
        assert len(orchestrator.ref_sequences) == 2
        assert len(orchestrator.ref_info) == 2
        
        # 验证参考序列内容
        assert 0 in orchestrator.ref_sequences
        assert 1 in orchestrator.ref_sequences
        assert orchestrator.ref_sequences[0].upper() == "ACGTACGTACGTACGTACGT"
        assert orchestrator.ref_sequences[1].upper() == "TGCATGCATGCATGCATGCA"
        
        # 验证参考信息
        assert orchestrator.ref_info[0][0] == "ref1"
        assert orchestrator.ref_info[0][1] == 20
        assert orchestrator.ref_info[1][0] == "ref2"
        assert orchestrator.ref_info[1][1] == 20
    
    @patch('sv_aligner.orchestrator.Seeder')
    @patch('sv_aligner.orchestrator.Chainer')
    @patch('sv_aligner.orchestrator.AlignerSW')
    @patch('sv_aligner.orchestrator.Refiner')
    @patch('sv_aligner.orchestrator.FinalSelector')
    def test_align_query(self, mock_selector, mock_refiner, mock_aligner, mock_chainer, mock_seeder):
        """测试查询序列比对流程"""
        # 配置mock对象，以返回我们期望的结果
        mock_seeder_instance = mock_seeder.return_value
        mock_seeder_instance.find_seeds.return_value = ["seed1", "seed2"]
        
        mock_chainer_instance = mock_chainer.return_value
        mock_chainer_instance.chain_seeds.return_value = ["chain1", "chain2"]
        
        mock_refiner_instance = mock_refiner.return_value
        mock_refiner_instance.refine_chains.return_value = [
            AlignmentSegment(q_name="query1", q_len=12, q_st=0, q_en=8, r_name="ref1", r_len=20, r_st=0, r_en=8, strand="+", score=16),
            AlignmentSegment(q_name="query1", q_len=12, q_st=8, q_en=12, r_name="ref1", r_len=20, r_st=8, r_en=12, strand="+", score=8)
        ]
        
        mock_selector_instance = mock_selector.return_value
        mock_selector_instance.select_alignments.return_value = [
            AlignmentSegment(q_name="query1", q_len=12, q_st=0, q_en=8, r_name="ref1", r_len=20, r_st=0, r_en=8, strand="+", score=16),
            AlignmentSegment(q_name="query1", q_len=12, q_st=8, q_en=12, r_name="ref1", r_len=20, r_st=8, r_en=12, strand="+", score=8)
        ]
        
        # 创建orchestrator实例并加载参考序列
        orchestrator = Orchestrator(self.params)
        orchestrator.ref_index = {}  # 假的索引，只为了不触发错误
        orchestrator.ref_info = {0: ("ref1", 20), 1: ("ref2", 20)}
        orchestrator.ref_sequences = {0: "ACGTACGTACGTACGTACGT", 1: "TGCATGCATGCATGCATGCA"}
        
        # 进行比对
        results = orchestrator.align_query(self.query_fasta)
        
        # 验证所有组件都被调用了
        assert mock_seeder.call_count == 1
        assert mock_chainer.call_count == 1
        assert mock_aligner.call_count == 1
        assert mock_refiner.call_count == 1
        assert mock_selector.call_count == 1
        
        # 验证返回了正确的结果
        assert len(results) == 2
        assert results[0].q_name == "query1"
        assert results[0].q_st == 0
        assert results[0].q_en == 8
        assert results[0].score == 16
        
        assert results[1].q_name == "query1"
        assert results[1].q_st == 8
        assert results[1].q_en == 12
        assert results[1].score == 8
    
    def test_align_query_real(self):
        """集成测试：使用真实组件进行完整的比对流程"""
        # 创建一个非常简单的查询和参考序列，以确保测试可以快速完成
        simple_ref_fasta = os.path.join(self.temp_dir.name, "simple_ref.fa")
        simple_query_fasta = os.path.join(self.temp_dir.name, "simple_query.fa")
        
        # 写入简单参考序列 - 只是一小段碱基
        with open(simple_ref_fasta, "w") as f:
            f.write(">simple_ref\n")
            f.write("ACGTACGTACGT\n")
        
        # 写入简单查询序列 - 完全匹配参考序列
        with open(simple_query_fasta, "w") as f:
            f.write(">simple_query\n")
            f.write("ACGTACGT\n")
        
        # 使用简化参数
        simple_params = dict(self.params)
        simple_params['k'] = 3  # 更小的k值让我们更容易找到种子
        simple_params['w'] = 1
        simple_params['min_segment_length'] = 4
        simple_params['min_alignment_score'] = 4
        simple_params['min_chain_score'] = 4
        
        # 创建orchestrator实例
        orchestrator = Orchestrator(simple_params)
        
        # 加载参考序列和索引
        orchestrator.load_reference_and_index(simple_ref_fasta)
        
        # 进行比对
        try:
            results = orchestrator.align_query(simple_query_fasta)
            
            # 验证找到了至少一个比对结果
            assert len(results) > 0
            
            # 验证第一个比对结果
            segment = results[0]
            assert segment.q_name == "simple_query"
            assert segment.q_len == 8
            assert segment.r_name == "simple_ref"
            assert segment.r_len == 12
            assert segment.strand == "+"
            assert segment.score > 0
            
        except Exception as e:
            # 如果发生异常，打印详细信息以帮助调试，但使测试失败
            print(f"实际比对测试失败：{str(e)}")
            raise
