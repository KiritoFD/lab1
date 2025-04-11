# SV-Aware 序列比对器

一个专为处理结构变异(SV)设计的高效DNA序列比对工具。

## 功能特点

该工具能够将查询DNA序列与参考DNA序列进行比对，特别关注以下类型的结构变异：

- 大型缺失 (Deletions)
- 大型插入 (Insertions)
- 倒置 (Inversions)
- 易位 (Translocations)
- 拷贝数变异/重复 (Duplications)

## 安装

### 前提条件

- Python 3.6 及以上版本
- parasail 库（用于 Smith-Waterman 比对）

### 安装步骤

```bash
# 克隆仓库
git clone https://github.com/yourusername/sv-aligner.git
cd sv-aligner

# 安装依赖项
pip install -r requirements.txt
```

## 使用方法

基本用法:

```bash
python run_sv_aligner.py reference.fa query.fa -o output.tsv
```

### 命令行参数

```
usage: run_sv_aligner.py [-h] [-o OUTPUT] [-k K] [-w W] [--index INDEX]
                         [--match MATCH_SCORE] [--mismatch MISMATCH_PENALTY]
                         [--gap-open GAP_OPEN_PENALTY]
                         [--gap-extend GAP_EXTEND_PENALTY]
                         [--min-len MIN_SEGMENT_LENGTH]
                         [--min-score MIN_ALIGNMENT_SCORE]
                         [--min-chain-score MIN_CHAIN_SCORE]
                         [--max-gap MAX_GAP]
                         [--max-dist-diff MAX_DIST_DIFF]
                         [--gap-penalty-factor GAP_PENALTY_FACTOR]
                         [--dist-diff-penalty-factor DIST_DIFF_PENALTY_FACTOR]
                         [-v]
                         reference query

位置参数:
  reference              参考基因组FASTA文件路径
  query                  查询序列FASTA文件路径

可选参数:
  -h, --help             显示帮助信息并退出
  -o OUTPUT, --output OUTPUT
                         输出TSV文件路径（默认：标准输出）
  -k K                   最小化子的k-mer大小（默认: 19）
  -w W                   最小化子窗口大小（默认: 10）
  --index INDEX          预构建索引的路径（如未提供或未找到，将构建索引）
  --match MATCH_SCORE    比对匹配得分（默认: 2）
  --mismatch MISMATCH_PENALTY
                         比对不匹配惩罚（默认: -3）
  --gap-open GAP_OPEN_PENALTY
                         比对缺口开启惩罚（默认: -5）
  --gap-extend GAP_EXTEND_PENALTY
                         比对缺口延伸惩罚（默认: -2）
  --min-len MIN_SEGMENT_LENGTH
                         要报告的比对段的最小长度（默认: 50）
  --min-score MIN_ALIGNMENT_SCORE
                         要报告的最小比对得分（默认: 40）
  --min-chain-score MIN_CHAIN_SCORE
                         考虑链接进行细化的最小链接得分（默认: 40）
  --max-gap MAX_GAP      链接期间允许的最大间隔大小（默认: 10000）
  --max-dist-diff MAX_DIST_DIFF
                         链接期间允许的查询和参考距离最大差异（默认: 500）
  --gap-penalty-factor GAP_PENALTY_FACTOR
                         链接期间间隔惩罚计算因子（默认: 0.01）
  --dist-diff-penalty-factor DIST_DIFF_PENALTY_FACTOR
                         链接期间距离差异惩罚因子（默认: 0.05）
  -v, --version          显示程序版本号并退出
```

## 输出格式

输出为制表符分隔值(TSV)格式，包含以下列：

1. `q_name`: 查询序列名称
2. `q_len`: 查询序列长度
3. `q_st`: 查询开始位置（0-based，包含）
4. `q_en`: 查询结束位置（0-based，不包含）
5. `r_name`: 参考序列名称
6. `r_len`: 参考序列长度
7. `r_st`: 参考开始位置（0-based，包含）
8. `r_en`: 参考结束位置（0-based，不包含）
9. `strand`: 链方向（'+'或'-'）
10. `score`: 比对得分
11. `edit_distance`: 编辑距离
12. `cigar`: CIGAR字符串（描述比对细节）

## 结构变异识别

可以通过分析连续的比对段之间的关系来识别结构变异：

- **缺失 (Deletion)**: 相邻行 i 和 j 有 `q_enᵢ ≈ q_stⱼ` 但 `r_stⱼ - r_enᵢ >> 0`
- **插入 (Insertion)**: 相邻行 i 和 j 有 `q_stⱼ - q_enᵢ >> 0` 但 `r_enᵢ ≈ r_stⱼ`
- **倒置 (Inversion)**: 行 i 的 `strand = '-'`
- **易位 (Translocation)**: 相邻行 i 和 j 有 `q_enᵢ ≈ q_stⱼ` 但参考位置跳跃很大或在不同染色体上
- **重复 (Duplication)**: 两个不同的行 i 和 j 映射到相同或大部分重叠的参考区域

## 算法概述

该比对器使用基于锚点(anchor-based)的策略：

1. **索引构建**: 使用最小化子(minimizers)为参考序列创建索引
2. **种子生成**: 通过查找查询序列与参考序列之间的共享最小化子来识别潜在的锚点
3. **链接**: 将共线性的种子链接成连贯的块
4. **精细化**: 使用Smith-Waterman算法对链接的区域进行局部比对
5. **选择**: 选择非重叠的高得分比对段作为最终结果

## 示例

```bash
# 比对查询序列到参考基因组并输出到文件
python run_sv_aligner.py reference.fasta query.fasta -o alignments.tsv

# 使用更严格的参数设置
python run_sv_aligner.py reference.fasta query.fasta --min-score 60 --min-len 100 -o strict_alignments.tsv

# 使用预构建的索引（加速重复运行）
python run_sv_aligner.py reference.fasta query.fasta --index ref_index.pkl -o alignments.tsv
```

## 项目结构

```
sv_aligner/
├── __init__.py          # 版本定义
├── data_types.py        # 核心数据类型定义
├── utils.py             # 工具函数
├── fasta_parser.py      # FASTA文件解析
├── indexer.py           # 参考基因组索引构建
├── seeder.py            # 种子生成
├── chainer.py           # 种子链接
├── aligner_sw.py        # Smith-Waterman实现
├── refiner.py           # 链精细化
├── selector.py          # 最终非重叠段选择
├── orchestrator.py      # 工作流协调
└── main.py              # 命令行接口
```

## 许可证

[在此添加许可信息]
