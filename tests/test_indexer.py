import os
import pytest
import tempfile
import pickle
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.indexer import IndexBuilder
from sv_aligner.data_types import ReferenceLocation

class TestIndexBuilder:
    def setup_method(self):
        """在每个测试方法前设置环境"""
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # 创建一个测试用的FASTA文件
        self.test_fasta = os.path.join(self.temp_dir.name, "test_ref.fa")
        with open(self.test_fasta, "w") as f:
            f.write(">ref1\n")
            f.write("ACGTACGTACGT\n")
            f.write(">ref2\n")
            f.write("TGCATGCATGCA\n")
        
        # 测试用的索引文件
        self.test_index = os.path.join(self.temp_dir.name, "test_index.pkl")
    
    def teardown_method(self):
        """在每个测试方法后清理环境"""
        self.temp_dir.cleanup()
    
    def test_build_index(self):
        """测试索引构建功能"""
        builder = IndexBuilder()
        ref_index, ref_info = builder.build(self.test_fasta, k=4, w=1)
        
        # 检查索引是否非空
        assert len(ref_index) > 0
        
        # 检查参考信息是否正确
        assert len(ref_info) == 2
        assert ref_info[0][0] == "ref1"
        assert ref_info[0][1] == 12
        assert ref_info[1][0] == "ref2"
        assert ref_info[1][1] == 12
        
        # 检查所有reference ID是否正确
        all_ref_ids = set()
        for locations in ref_index.values():
            for loc in locations:
                all_ref_ids.add(loc.ref_id)
        
        assert all_ref_ids == {0, 1}  # 应该包含两个参考序列的ID
    
    def test_save_load_index(self):
        """测试索引保存和加载功能"""
        builder = IndexBuilder()
        ref_index_original, ref_info_original = builder.build(self.test_fasta, k=4, w=1)
        
        # 保存索引
        builder.save_index(self.test_index)
        
        # 检查索引文件是否被创建
        assert os.path.exists(self.test_index)
        
        # 从文件加载索引
        loaded_index, loaded_info = IndexBuilder.load_index(self.test_index)
        
        # 验证加载的索引与原始索引相同
        assert len(loaded_index) == len(ref_index_original)
        assert loaded_info == ref_info_original
    
    def test_empty_sequence(self):
        """测试处理空序列"""
        empty_fasta = os.path.join(self.temp_dir.name, "empty.fa")
        with open(empty_fasta, "w") as f:
            f.write(">empty1\n\n")
            f.write(">empty2\n\n")
        
        builder = IndexBuilder()
        ref_index, ref_info = builder.build(empty_fasta, k=4, w=1)
        
        # 空序列应该被跳过，但信息应该被记录
        assert len(ref_index) == 0
        assert len(ref_info) == 0
    
    def test_load_nonexistent_index(self):
        """测试加载不存在的索引文件"""
        non_existent = os.path.join(self.temp_dir.name, "non_existent.pkl")
        with pytest.raises(FileNotFoundError):
            IndexBuilder.load_index(non_existent)
