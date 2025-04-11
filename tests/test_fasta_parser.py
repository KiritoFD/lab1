import os
import pytest
import tempfile
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.fasta_parser import read_fasta

class TestFastaParser:
    def setup_method(self):
        """在每个测试方法前创建临时FASTA文件"""
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # 创建一个简单的FASTA文件
        self.simple_fasta = os.path.join(self.temp_dir.name, "simple.fa")
        with open(self.simple_fasta, "w") as f:
            f.write(">seq1\n")
            f.write("ACGTACGT\n")
            f.write(">seq2\n")
            f.write("TGCATGCA\n")
        
        # 创建一个多行序列的FASTA文件
        self.multiline_fasta = os.path.join(self.temp_dir.name, "multiline.fa")
        with open(self.multiline_fasta, "w") as f:
            f.write(">seq1 description here\n")
            f.write("ACGT\n")
            f.write("ACGT\n")
            f.write(">seq2 another description\n")
            f.write("TGCA\n")
            f.write("TGCA\n")
        
        # 创建一个空FASTA文件
        self.empty_fasta = os.path.join(self.temp_dir.name, "empty.fa")
        open(self.empty_fasta, "w").close()
        
        # 创建一个无效FASTA文件（没有header）
        self.invalid_fasta = os.path.join(self.temp_dir.name, "invalid.fa")
        with open(self.invalid_fasta, "w") as f:
            f.write("ACGTACGT\n")
    
    def teardown_method(self):
        """在每个测试方法后删除临时文件"""
        self.temp_dir.cleanup()
    
    def test_read_simple_fasta(self):
        """测试读取简单FASTA文件"""
        sequences = list(read_fasta(self.simple_fasta))
        assert len(sequences) == 2
        assert sequences[0] == ("seq1", "ACGTACGT")
        assert sequences[1] == ("seq2", "TGCATGCA")
    
    def test_read_multiline_fasta(self):
        """测试读取多行序列的FASTA文件"""
        sequences = list(read_fasta(self.multiline_fasta))
        assert len(sequences) == 2
        # 应该只取header的第一部分作为名称
        assert sequences[0][0] == "seq1"
        assert sequences[0][1] == "ACGTACGT"
        assert sequences[1][0] == "seq2"
        assert sequences[1][1] == "TGCATGCA"
    
    def test_read_empty_fasta(self):
        """测试读取空FASTA文件"""
        sequences = list(read_fasta(self.empty_fasta))
        assert len(sequences) == 0
    
    def test_read_invalid_fasta(self):
        """测试读取无效FASTA文件"""
        sequences = list(read_fasta(self.invalid_fasta))
        assert len(sequences) == 0  # 应该不返回任何序列，因为没有header
    
    def test_file_not_found(self):
        """测试读取不存在的文件"""
        non_existent = os.path.join(self.temp_dir.name, "non_existent.fa")
        with pytest.raises(FileNotFoundError):
            list(read_fasta(non_existent))
