import sys
from typing import Tuple, Optional, Dict, Any

# Try to import parasail, but have a fallback for testing
try:
    import parasail
    PARASAIL_AVAILABLE = True
except ImportError:
    PARASAIL_AVAILABLE = False
    print("Warning: parasail library not found. Using native Python Smith-Waterman implementation (slow).", file=sys.stderr)
    
    # Native implementation for testing when parasail isn't available
    def native_sw(query, ref, match_score, mismatch_penalty, gap_open, gap_extend):
        """Very basic Smith-Waterman implementation for testing purposes."""
        m, n = len(query), len(ref)
        score_matrix = [[0 for _ in range(n+1)] for _ in range(m+1)]
        backtrack = [[0 for _ in range(n+1)] for _ in range(m+1)]  # 0: stop, 1: diag, 2: up, 3: left
        max_score, max_i, max_j = 0, 0, 0
        
        for i in range(1, m+1):
            for j in range(1, n+1):
                match = score_matrix[i-1][j-1] + (match_score if query[i-1] == ref[j-1] else mismatch_penalty)
                delete = score_matrix[i-1][j] + gap_open if backtrack[i-1][j] != 2 else score_matrix[i-1][j] + gap_extend
                insert = score_matrix[i][j-1] + gap_open if backtrack[i][j-1] != 3 else score_matrix[i][j-1] + gap_extend
                
                score_matrix[i][j] = max(0, match, delete, insert)
                
                if score_matrix[i][j] == 0:
                    backtrack[i][j] = 0
                elif score_matrix[i][j] == match:
                    backtrack[i][j] = 1
                elif score_matrix[i][j] == delete:
                    backtrack[i][j] = 2
                else:
                    backtrack[i][j] = 3
                    
                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_i, max_j = i, j
                    
        # Trace back to find start positions and build CIGAR
        if max_score == 0:
            return 0, "0M", 0, 0, 0, 0
            
        i, j = max_i, max_j
        q_start, r_start = i, j
        cigar_ops = []
        current_op = None
        count = 0
        
        while i > 0 and j > 0 and backtrack[i][j] != 0:
            if backtrack[i][j] == 1:  # diagonal - match/mismatch
                op = 'M'
                i -= 1
                j -= 1
            elif backtrack[i][j] == 2:  # up - deletion in reference
                op = 'I'
                i -= 1
            else:  # left - insertion in reference
                op = 'D'
                j -= 1
                
            if op != current_op:
                if current_op:
                    cigar_ops.append(f"{count}{current_op}")
                current_op = op
                count = 1
            else:
                count += 1
                
        if count > 0:
            cigar_ops.append(f"{count}{current_op}")
            
        cigar_string = "".join(reversed(cigar_ops))
        q_start = i
        r_start = j
        
        return max_score, cigar_string, q_start, max_i, r_start, max_j

class AlignerSW:
    def __init__(self, params: Dict[str, Any]):
        self.match_score = params.get('match_score', 2)
        self.mismatch_penalty = params.get('mismatch_penalty', -3)
        self.gap_open_penalty = params.get('gap_open_penalty', -5)
        self.gap_extend_penalty = params.get('gap_extend_penalty', -2)
        
        if PARASAIL_AVAILABLE:
            # Create substitution matrix (simple for DNA)
            self.scoring_matrix = parasail.matrix_create("ACGT", self.match_score, -self.mismatch_penalty)

    def align(self, query_sub: str, ref_sub: str) -> Optional[Tuple[int, str, int, int, int, int]]:
        """
        Performs Smith-Waterman alignment using Parasail or native implementation.

        Returns:
            Tuple (score, cigar_string, query_begin, query_end, ref_begin, ref_end) or None
            Coordinates are 0-based relative to the input subsequences.
        """
        if not query_sub or not ref_sub:
            return None

        try:
            if PARASAIL_AVAILABLE:
                # Use parasail.sg for Smith-Waterman Gotoh (affine gap)
                result = parasail.sg_stats_striped_32(
                    query_sub,
                    ref_sub,
                    abs(self.gap_open_penalty),  # Parasail uses positive penalties
                    abs(self.gap_extend_penalty),
                    self.scoring_matrix
                )

                score = result.score
                # End positions are inclusive in parasail, convert to exclusive for q_end, r_end
                q_end = result.end_query + 1
                r_end = result.end_ref + 1

                # Get CIGAR string
                cigar = result.cigar.decode

                # Estimate start positions based on CIGAR
                q_len_aligned = result.cigar.query_length
                r_len_aligned = result.cigar.ref_length
                q_st = max(0, q_end - q_len_aligned)
                r_st = max(0, r_end - r_len_aligned)

                return score, cigar, q_st, q_end, r_st, r_end
            else:
                # Use native implementation if Parasail is not available
                return native_sw(
                    query_sub, 
                    ref_sub, 
                    self.match_score,
                    self.mismatch_penalty,
                    self.gap_open_penalty,
                    self.gap_extend_penalty
                )

        except Exception as e:
            print(f"Error during alignment: {e}", file=sys.stderr)
            return None
