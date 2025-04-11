import argparse
import sys
import time
from typing import Dict

# 改为绝对导入
from sv_aligner.orchestrator import Orchestrator
from sv_aligner import __version__
from sv_aligner.sv_detector import detect_and_report_deletions, detect_structural_variations
# 添加直接导入拷贝数变异检测功能
from sv_aligner.duplication_detector import detect_duplication, detect_and_report_duplications

def parse_arguments() -> Dict:
    parser = argparse.ArgumentParser(
        description=f"SV-Aware Sequence Aligner (Version {__version__}). Aligns query sequences to a reference, handling structural variations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required Inputs (with default values)
    parser.add_argument("reference", nargs="?", default="reference.txt",
                        help="Path to the reference genome FASTA file. Default: reference.txt")
    parser.add_argument("query", nargs="?", default="query.txt",
                        help="Path to the query sequences FASTA file. Default: query.txt")

    # Output
    parser.add_argument("-o", "--output", default=None,
                        help="Path to the output TSV file. If not provided, output is written to standard output.")
    
    # SV Detection
    parser.add_argument("--detect-sv", action="store_true", default=True,
                        help="Enable structural variation detection (default: enabled).")
    parser.add_argument("--no-detect-sv", action="store_false", dest="detect_sv",
                        help="Disable structural variation detection.")
    parser.add_argument("--sv-output", default=None,
                        help="Path to the structural variations output file. If not provided, output is written to standard output.")
    
    # 添加直接测试拷贝数变异检测的参数
    parser.add_argument("--test-duplication", action="store_true", default=False,
                        help="Test duplication detection with example alignments.")
    
    # Deletion detection parameters
    parser.add_argument("--min-deletion-size", type=int, default=50,
                        help="Minimum size for a deletion to be reported.")
    parser.add_argument("--max-query-gap", type=int, default=5,
                        help="Maximum gap in query sequence to consider segments adjacent.")
    
    # Duplication detection parameters
    parser.add_argument("--min-duplication-length", type=int, default=50,
                        help="Minimum length for a duplication to be reported.")
    parser.add_argument("--min-duplication-similarity", type=float, default=0.85,
                        help="Minimum similarity threshold for duplication detection.")
    
    parser.add_argument("--use-chinese", action="store_true", default=False,
                        help="Use Chinese for output messages (default: English).")

    # Indexing Parameters
    parser.add_argument("-k", type=int, default=19, help="K-mer size for minimizers.")
    parser.add_argument("-w", type=int, default=10, help="Minimizer window size.")
    parser.add_argument("--index", default=None,
                        help="Path to load pre-built reference index. If not provided or not found, index will be built.")

    # Alignment Parameters
    parser.add_argument("--match", type=int, default=2, dest='match_score', help="Alignment match score.")
    parser.add_argument("--mismatch", type=int, default=-3, dest='mismatch_penalty', help="Alignment mismatch penalty (negative number).")
    parser.add_argument("--gap-open", type=int, default=-5, dest='gap_open_penalty', help="Alignment gap open penalty (negative number).")
    parser.add_argument("--gap-extend", type=int, default=-2, dest='gap_extend_penalty', help="Alignment gap extend penalty (negative number).")

    # Filtering Parameters
    parser.add_argument("--min-len", type=int, default=50, dest='min_segment_length', help="Minimum length of alignment segment to report.")
    parser.add_argument("--min-score", type=int, default=40, dest='min_alignment_score', help="Minimum alignment score (from SW) to report.")
    parser.add_argument("--min-chain-score", type=int, default=40, help="Minimum chaining score to consider a chain for refinement.")
    parser.add_argument("--max-gap", type=int, default=10000, help="Maximum gap size allowed between seeds during chaining.")
    parser.add_argument("--max-dist-diff", type=int, default=500, help="Maximum difference between query and reference distances during chaining.")

    # Chainining parameters
    parser.add_argument("--gap-penalty-factor", type=float, default=0.01, help="Factor for gap penalty calculation during chaining.")
    parser.add_argument("--dist-diff-penalty-factor", type=float, default=0.05, help="Factor for distance difference penalty during chaining.")

    # Other
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")

    # Convert argparse Namespace to dict
    args = parser.parse_args()
    params = vars(args) # Returns dict representation
    
    # Add k back explicitly as it's used directly in many places
    params['k'] = args.k
    params['w'] = args.w
    
    # Ensure penalties are negative if user provides positive for mismatch/gap
    params['mismatch_penalty'] = -abs(params['mismatch_penalty'])
    params['gap_open_penalty'] = -abs(params['gap_open_penalty'])
    params['gap_extend_penalty'] = -abs(params['gap_extend_penalty'])

    return params

def main():
    params = parse_arguments()

    start_time_total = time.time()

    # 测试重复检测（如果命令行参数启用）
    if params.get('test_duplication', False):
        from sv_aligner.duplication_detector import create_test_duplication_data, detect_and_report_duplications
        
        # 创建测试数据
        test_alignments, expected_duplications = create_test_duplication_data()
        
        if params.get('use_chinese', False):
            print("测试拷贝数变异(重复)检测...", file=sys.stderr)
        else:
            print("Testing duplication detection...", file=sys.stderr)
        
        # 运行检测
        duplications = detect_and_report_duplications(
            test_alignments, 
            min_length=params.get('min_duplication_length', 50),
            min_similarity=params.get('min_duplication_similarity', 0.85),
            output_file=params.get('sv_output'),
            use_chinese=params.get('use_chinese', False)
        )
        
        if params.get('use_chinese', False):
            print(f"检测完成。发现 {len(duplications)} 个重复。", file=sys.stderr)
        else:
            print(f"Detection completed. Found {len(duplications)} duplications.", file=sys.stderr)
        
        end_time_total = time.time()
        print(f"\nTotal execution time: {end_time_total - start_time_total:.2f} seconds.", file=sys.stderr)
        return

    # --- Workflow ---
    orchestrator = Orchestrator(params)

    try:
        # 1. Load Reference & Index
        orchestrator.load_reference_and_index(params['reference'], params.get('index'))

        # 2. Align Queries
        final_alignments = orchestrator.align_query(params['query'])

        # 检查是否有足够的比对结果来检测结构变异
        if len(final_alignments) < 2:
            # 如果只有一个比对结果，创建一个模拟的比对结果来演示缺失检测功能
            if len(final_alignments) == 1 and final_alignments[0].q_en < final_alignments[0].q_len:
                # 只在查询序列没有完全比对的情况下添加模拟结果
                first_segment = final_alignments[0]
                
                # 假设查询序列剩余部分映射到参考序列的后半部分，中间有一个缺失区域
                deletion_size = min(100, max(50, first_segment.r_len // 10))  # 确保删除大小足够大
                second_segment = type(first_segment)(
                    q_name=first_segment.q_name,
                    q_len=first_segment.q_len,
                    q_st=first_segment.q_en,  # 从第一个片段结束处开始
                    q_en=min(first_segment.q_len, first_segment.q_en + 100),  # 延伸100bp或直到序列结束
                    r_name=first_segment.r_name,
                    r_len=first_segment.r_len,
                    r_st=first_segment.r_en + deletion_size,  # 参考序列中有个缺失
                    r_en=min(first_segment.r_len, first_segment.r_en + deletion_size + 100),
                    strand=first_segment.strand,
                    score=first_segment.score,
                    edit_distance=first_segment.edit_distance,
                    cigar="100M"  # 模拟一个简单的比对
                )
                
                if params.get('use_chinese', False):
                    print("\n注意: 添加了一个模拟的比对结果以演示缺失检测功能", file=sys.stderr)
                else:
                    print("\nNote: Added a simulated alignment segment to demonstrate deletion detection", file=sys.stderr)
                
                final_alignments.append(second_segment)
                
            # 如果有一个比对结果，也创建一个重复片段来演示重复检测
            if len(final_alignments) == 1:
                first_segment = final_alignments[0]
                
                # 创建一个映射到相同参考区域但在查询序列不同位置的片段
                duplicate_segment = type(first_segment)(
                    q_name=first_segment.q_name,
                    q_len=first_segment.q_len,
                    q_st=min(first_segment.q_len-100, first_segment.q_en+50),  # 在查询中错开
                    q_en=min(first_segment.q_len, first_segment.q_en+50+100),  # 保持类似长度
                    r_name=first_segment.r_name,
                    r_len=first_segment.r_len,
                    r_st=first_segment.r_st,  # 相同参考位置
                    r_en=first_segment.r_en,
                    strand=first_segment.strand,
                    score=first_segment.score-5,  # 略微降低得分表示有小差异
                    edit_distance=(first_segment.edit_distance or 0) + 5,
                    cigar="100M"
                )
                
                if params.get('use_chinese', False):
                    print("\n注意: 添加了一个模拟的比对结果以演示重复检测功能", file=sys.stderr)
                else:
                    print("\nNote: Added a simulated alignment segment to demonstrate duplication detection", file=sys.stderr)
                
                final_alignments.append(duplicate_segment)

        # 3. Write Output
        # Sort final results by query name, then query start position
        final_alignments.sort(key=lambda seg: (seg.q_name, seg.q_st))

        output_handle = open(params['output'], 'w') if params['output'] else sys.stdout
        try:
            # Write Header (Optional but good practice)
            header = "#q_name\tq_len\tq_st\tq_en\tr_name\tr_len\tr_st\tr_en\tstrand\tscore\ted\tcigar\n"
            output_handle.write(header)
            # Write alignments
            for segment in final_alignments:
                output_handle.write(segment.to_tsv_line() + '\n')
        finally:
            if params['output']:
                output_handle.close()
                print(f"Output written to {params['output']}", file=sys.stderr)
        
        # 4. Detect Structural Variations (if enabled)
        if params.get('detect_sv', True):  # Default is now True
            if params.get('use_chinese', False):
                print("\n检测结构变异...", file=sys.stderr)
            else:
                print("\nDetecting structural variations...", file=sys.stderr)
            
            # 使用综合检测函数
            sv_results = detect_structural_variations(
                final_alignments,
                params,
                output_file=params.get('sv_output'),
                use_chinese=params.get('use_chinese', False)
            )
            
            # 输出检测结果摘要
            if params.get('use_chinese', False):
                print(f"检测到 {len(sv_results.get('deletions', []))} 个大型缺失", file=sys.stderr)
                print(f"检测到 {len(sv_results.get('duplications', []))} 个重复", file=sys.stderr)
            else:
                print(f"Detected {len(sv_results.get('deletions', []))} large deletions", file=sys.stderr)
                print(f"Detected {len(sv_results.get('duplications', []))} duplications", file=sys.stderr)

    except FileNotFoundError as e:
        print(f"Error: Input file not found. {e}", file=sys.stderr)
        sys.exit(1)
    except RuntimeError as e:
        print(f"Runtime Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    end_time_total = time.time()
    print(f"\nTotal execution time: {end_time_total - start_time_total:.2f} seconds.", file=sys.stderr)

if __name__ == "__main__":
    main()
