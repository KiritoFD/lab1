#!/usr/bin/env python3
"""
测试脚本：使用ref.txt作为参考序列，que.txt作为查询序列进行测试
"""

import os
import sys
import subprocess
import tempfile
from pathlib import Path

def create_fasta_from_txt(txt_path, output_fasta_path, header_name):
    """将txt格式的序列文件转换为FASTA格式"""
    try:
        # 读取txt文件内容
        with open(txt_path, 'r') as f:
            sequence = f.read().strip()
        
        # 写入FASTA格式
        with open(output_fasta_path, 'w') as f:
            f.write(f">{header_name}\n")
            # 每行最多80个字符
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")
                
        print(f"已创建FASTA文件: {output_fasta_path}")
        return True
    except Exception as e:
        print(f"转换文件时出错: {str(e)}")
        return False

def run_test(ref_fasta, query_fasta, output_file):
    """运行SV-Aware序列比对器测试"""
    try:
        # 构建命令
        cmd = [
            sys.executable, 
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_sv_aligner.py"),
            ref_fasta,
            query_fasta,
            "-o", output_file,
            # 调整参数以适应测试数据
            "--min-len", "30",
            "--min-score", "20",
            "--min-chain-score", "20",
            "-k", "15",
            "-w", "5"
        ]
        
        # 执行命令
        print(f"执行命令: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # 检查执行结果
        if result.returncode == 0:
            print("比对器成功运行")
            print("\n--- 标准输出 ---")
            print(result.stdout)
        else:
            print(f"比对器运行失败，返回码: {result.returncode}")
            print("\n--- 标准错误 ---")
            print(result.stderr)
            return False
            
        # 检查输出文件
        if os.path.exists(output_file):
            print(f"\n输出文件 '{output_file}' 已创建")
            with open(output_file, 'r') as f:
                print("\n--- 比对结果 ---")
                print(f.read())
        else:
            print(f"错误: 输出文件 '{output_file}' 未创建")
            return False
            
        return True
    except Exception as e:
        print(f"测试执行时出错: {str(e)}")
        return False

def main():
    # 项目根目录
    project_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 输入文件路径
    ref_txt_path = os.path.join(project_dir, "python", "ref.txt")
    que_txt_path = os.path.join(project_dir, "python", "que.txt")
    
    # 检查输入文件是否存在
    if not os.path.exists(ref_txt_path):
        print(f"错误: 参考文件 '{ref_txt_path}' 不存在")
        return 1
        
    # 为测试创建查询文件（如果不存在）
    if not os.path.exists(que_txt_path):
        print(f"查询文件 '{que_txt_path}' 不存在，创建测试查询文件...")
        try:
            # 从ref.txt中提取子序列创建查询序列
            with open(ref_txt_path, 'r') as f:
                ref_seq = f.read().strip()
                
            # 创建包含结构变异的查询序列
            # 1. 简单匹配段
            query_seq = ref_seq[100:600]
            # 2. 倒置段
            inv_start, inv_end = 700, 900
            inv_seq = ref_seq[inv_start:inv_end]
            rev_comp = "".join({'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}.get(base.upper(), base) 
                             for base in reversed(inv_seq))
            query_seq += rev_comp
            # 3. 另一个匹配段（模拟易位）
            query_seq += ref_seq[1500:1800]
            
            # 写入que.txt文件
            with open(que_txt_path, 'w') as f:
                f.write(query_seq)
                
            print(f"已创建测试查询文件: {que_txt_path} (长度: {len(query_seq)})")
        except Exception as e:
            print(f"创建查询文件时出错: {str(e)}")
            return 1
    
    # 创建临时FASTA文件
    ref_fasta = os.path.join(project_dir, "ref.fa")
    query_fasta = os.path.join(project_dir, "query.fa")
    
    # 转换txt到FASTA格式
    if not create_fasta_from_txt(ref_txt_path, ref_fasta, "reference"):
        return 1
    if not create_fasta_from_txt(que_txt_path, query_fasta, "query"):
        return 1
    
    # 设置输出文件路径
    output_file = os.path.join(project_dir, "alignment_results.tsv")
    
    # 运行测试
    success = run_test(ref_fasta, query_fasta, output_file)
    
    # 清理临时文件（可选）
    # os.remove(ref_fasta)
    # os.remove(query_fasta)
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
