#!/usr/bin/env python3
"""
将纯序列文件转换为FASTA格式
"""
import os
import sys
import argparse

def convert_to_fasta(input_file, output_file=None, sequence_name=None):
    """
    将纯序列文件转换为FASTA格式
    
    Args:
        input_file: 输入文件路径
        output_file: 输出FASTA文件路径，如果不提供则使用输入文件名加上.fa后缀
        sequence_name: 序列名称，如果不提供则使用文件名
    """
    if not os.path.exists(input_file):
        print(f"错误: 输入文件不存在: {input_file}", file=sys.stderr)
        return False
        
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}.fa"
        
    if sequence_name is None:
        sequence_name = os.path.basename(os.path.splitext(input_file)[0])
    
    try:
        with open(input_file, 'r') as infile:
            sequence_data = []
            for line in infile:
                line = line.strip()
                if line:  # 跳过空行
                    sequence_data.append(line)
            
            if not sequence_data:
                print(f"错误: 输入文件为空: {input_file}", file=sys.stderr)
                return False
            
            sequence = "".join(sequence_data)
            
            # 检查是否为有效的DNA序列
            valid_chars = set("ACGTNacgtn")
            invalid_chars = set(sequence) - valid_chars
            if invalid_chars:
                print(f"警告: 序列包含非标准DNA字符: {', '.join(invalid_chars)}", file=sys.stderr)
                print(f"将继续处理，但非标准字符可能导致问题", file=sys.stderr)
        
        with open(output_file, 'w') as outfile:
            outfile.write(f">{sequence_name}\n")
            
            # 每行最多80个字符
            for i in range(0, len(sequence), 80):
                outfile.write(sequence[i:i+80] + "\n")
        
        print(f"成功转换 {input_file} 为FASTA格式")
        print(f"输出文件: {output_file}")
        print(f"序列名称: {sequence_name}")
        print(f"序列长度: {len(sequence)} bp")
        return True
            
    except Exception as e:
        print(f"错误: 转换过程中出现异常: {e}", file=sys.stderr)
        return False

def main():
    parser = argparse.ArgumentParser(description="将纯序列文件转换为FASTA格式")
    parser.add_argument("input_file", help="输入序列文件路径")
    parser.add_argument("-o", "--output", help="输出FASTA文件路径")
    parser.add_argument("-n", "--name", help="序列名称")
    
    args = parser.parse_args()
    
    success = convert_to_fasta(args.input_file, args.output, args.name)
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
