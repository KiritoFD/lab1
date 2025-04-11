import sys
import os
from typing import Tuple, Iterator

def read_fasta(filepath: str) -> Iterator[Tuple[str, str]]:
    """
    读取FASTA文件或纯序列文件，返回(header, sequence)元组的迭代器。
    处理多行序列。
    """
    # 检查文件是否存在
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"序列文件不存在: {filepath}")
        
    header = None
    sequence = []
    try:
        with open(filepath, 'r') as infile:
            # 检查文件是否为空
            is_empty = True
            is_fasta = False
            
            for line_num, line in enumerate(infile):
                is_empty = False
                line = line.strip()
                if not line:
                    continue
                
                # 如果是第一行且以'>'开头，则认为是FASTA格式
                if line_num == 0 and line.startswith('>'):
                    is_fasta = True
                    header = line[1:].split()[0] if ' ' in line[1:] else line[1:]
                    continue
                
                # 对于FASTA格式的处理
                if is_fasta:
                    if line.startswith('>'):
                        if header is not None:
                            yield header, "".join(sequence)
                        header = line[1:].split()[0] if ' ' in line[1:] else line[1:]
                        sequence = []
                    elif header is not None:  # 确保我们正在处理序列行
                        sequence.append(line)
                # 对于纯序列文件的处理
                else:
                    # 检查是否为有效的DNA序列字符
                    if all(c.upper() in "ACGTN" for c in line):
                        sequence.append(line)
                    else:
                        print(f"警告: 文件 {filepath} 第 {line_num+1} 行包含非DNA序列字符", file=sys.stderr)
            
            # 输出最后一个序列
            if is_fasta:
                if header is not None:
                    yield header, "".join(sequence)
            else:
                if sequence:
                    # 对于纯序列文件，使用文件名作为序列名称
                    filename = os.path.basename(filepath)
                    yield filename.split('.')[0], "".join(sequence)
            
            # 检查是否读取到了有效的序列
            if is_empty or (is_fasta and header is None) or (not is_fasta and not sequence):
                print(f"警告: 文件 {filepath} 似乎不包含有效的序列数据", file=sys.stderr)
                
    except UnicodeDecodeError:
        raise RuntimeError(f"文件 {filepath} 不是有效的文本文件。请确保它是ASCII或UTF-8编码的序列文件。")
    except Exception as e:
        raise RuntimeError(f"读取序列文件 {filepath} 时出错: {e}")
