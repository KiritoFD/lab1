#!/usr/bin/env python3
"""
针对纯序列文件的SV-Aligner运行脚本
"""
import os
import sys
import argparse
import subprocess
import tempfile

def convert_to_fasta_temp(seq_file, seq_name=None):
    """将序列文件转换为临时FASTA文件"""
    if not os.path.exists(seq_file):
        print(f"错误: 序列文件不存在: {seq_file}", file=sys.stderr)
        return None
    
    if seq_name is None:
        seq_name = os.path.basename(os.path.splitext(seq_file)[0])
    
    try:
        # 创建临时文件
        temp_fd, temp_path = tempfile.mkstemp(suffix=".fa")
        os.close(temp_fd)
        
        with open(seq_file, 'r') as infile:
            sequence_data = []
            for line in infile:
                line = line.strip()
                if line:  # 跳过空行
                    sequence_data.append(line)
            
            if not sequence_data:
                print(f"错误: 序列文件为空: {seq_file}", file=sys.stderr)
                os.unlink(temp_path)
                return None
            
            sequence = "".join(sequence_data)
        
        with open(temp_path, 'w') as outfile:
            outfile.write(f">{seq_name}\n")
            
            # 每行最多80个字符
            for i in range(0, len(sequence), 80):
                outfile.write(sequence[i:i+80] + "\n")
        
        print(f"已将 {seq_file} 转换为临时FASTA文件")
        print(f"序列名称: {seq_name}")
        print(f"序列长度: {len(sequence)} bp")
        return temp_path
            
    except Exception as e:
        print(f"错误: 转换过程中出现异常: {e}", file=sys.stderr)
        if 'temp_path' in locals():
            os.unlink(temp_path)
        return None

def main():
    parser = argparse.ArgumentParser(description="使用纯序列文件运行SV-Aligner")
    parser.add_argument("reference", nargs="?", default="ref.txt",
                        help="参考序列文件路径 (默认: ref.txt)")
    parser.add_argument("query", nargs="?", default="que.txt",
                        help="查询序列文件路径 (默认: que.txt)")
    parser.add_argument("-o", "--output", default=None,
                        help="输出文件路径")
    parser.add_argument("-r", "--ref-name", default="reference",
                        help="参考序列名称 (默认: reference)")
    parser.add_argument("-q", "--query-name", default="query",
                        help="查询序列名称 (默认: query)")
    parser.add_argument("--keep-temp", action="store_true",
                        help="保留临时FASTA文件")
    
    # 添加其他SV-Aligner参数
    parser.add_argument("-k", type=int, default=19,
                        help="K-mer大小 (默认: 19)")
    parser.add_argument("-w", type=int, default=10,
                        help="最小化子窗口大小 (默认: 10)")
    
    args = parser.parse_args()
    
    # 转换参考序列文件
    ref_fasta = convert_to_fasta_temp(args.reference, args.ref_name)
    if not ref_fasta:
        return 1
    
    # 转换查询序列文件
    query_fasta = convert_to_fasta_temp(args.query, args.query_name)
    if not query_fasta:
        os.unlink(ref_fasta)
        return 1
    
    try:
        # 构建SV-Aligner命令
        cmd = [sys.executable, "-m", "sv_aligner.main"]
        cmd.extend([ref_fasta, query_fasta])
        
        if args.output:
            cmd.extend(["-o", args.output])
        
        # 添加其他参数
        cmd.extend(["-k", str(args.k)])
        cmd.extend(["-w", str(args.w)])
        
        # 运行SV-Aligner
        print(f"运行命令: {' '.join(cmd)}", file=sys.stderr)
        result = subprocess.run(cmd)
        
        return result.returncode
    
    finally:
        # 清理临时文件
        if not args.keep_temp:
            print("清理临时文件...", file=sys.stderr)
            os.unlink(ref_fasta)
            os.unlink(query_fasta)
        else:
            print(f"保留临时文件: {ref_fasta}, {query_fasta}", file=sys.stderr)

if __name__ == "__main__":
    sys.exit(main())
