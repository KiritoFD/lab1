import sys
from typing import List, Dict, Tuple, Any

from .data_types import Seed, AlignmentSegment
from .utils import reverse_complement
from .aligner_sw import AlignerSW

class Refiner:
    def __init__(self, aligner: AlignerSW, ref_sequences: Dict[int, str], ref_info: Dict[int, Tuple[str, int]], params: Dict[str, Any]):
        self.aligner = aligner
        self.ref_sequences = ref_sequences # Need actual sequences {ref_id: seq_str}
        self.ref_info = ref_info # {ref_id: (name, len)}
        self.params = params

    def refine_chains(self, chains: List[List[Seed]], query_seq: str, q_name: str) -> List[AlignmentSegment]:
        """Refines seed chains into detailed alignment segments."""
        refined_segments = []
        q_len = len(query_seq)
        min_segment_len = self.params.get('min_segment_length', 50)
        min_alignment_score = self.params.get('min_alignment_score', 30)

        for chain_idx, chain in enumerate(chains):
            if not chain: continue

            # --- Estimate Region ---
            # Use first and last seed positions as rough guides
            first_seed = chain[0]
            last_seed = chain[-1]
            chain_strand = first_seed.strand # Strand is consistent within a chain

            # Add padding (e.g., k-mer length or a fixed amount)
            padding = self.params.get('k', 15) * 2 # Example padding

            q_start_est = max(0, first_seed.q_pos - padding)
            q_end_est = min(q_len, last_seed.q_pos + self.params.get('k', 15) + padding)

            ref_id = first_seed.ref_id
            ref_name, ref_len = self.ref_info[ref_id]
            r_start_est = max(0, first_seed.r_pos - padding)
            r_end_est = min(ref_len, last_seed.r_pos + self.params.get('k', 15) + padding)

            if q_start_est >= q_end_est or r_start_est >= r_end_est:
                continue

            # --- Extract Sequences ---
            query_sub = query_seq[q_start_est:q_end_est]
            ref_seq_full = self.ref_sequences.get(ref_id)
            if not ref_seq_full:
                 print(f"Error: Could not find reference sequence for ref_id {ref_id}", file=sys.stderr)
                 continue # Should not happen if index built correctly

            ref_sub_oriented = ref_seq_full[r_start_est:r_end_est]
            if not chain_strand: # If strand is False ('-'), use reverse complement
                ref_sub_oriented = reverse_complement(ref_sub_oriented)

            if not query_sub or not ref_sub_oriented:
                continue

            # --- Align ---
            alignment_result = self.aligner.align(query_sub, ref_sub_oriented)

            if alignment_result:
                score, cigar, q_sub_st, q_sub_en, r_sub_st, r_sub_en = alignment_result

                # --- Convert back to original coordinates ---
                final_q_st = q_start_est + q_sub_st
                final_q_en = q_start_est + q_sub_en

                # Reference coordinates need care depending on strand
                if chain_strand: # '+' strand
                    final_r_st = r_start_est + r_sub_st
                    final_r_en = r_start_est + r_sub_en
                else: # '-' strand
                    # Coordinates from alignment are relative to the RC subsequence
                    # Map back to original reference forward strand coordinates
                    original_len_r_sub = r_end_est - r_start_est
                    final_r_st = r_end_est - r_sub_en # Start on forward strand = end of original - end on RC
                    final_r_en = r_end_est - r_sub_st # End on forward strand = end of original - start on RC

                # --- Filter ---
                aln_len = final_q_en - final_q_st # Use query length as representative
                if aln_len >= min_segment_len and score >= min_alignment_score:
                    # Calculate edit distance from CIGAR
                    edit_distance = self._calculate_edit_distance_from_cigar(cigar)
                    
                    segment = AlignmentSegment(
                        q_name=q_name, q_len=q_len, q_st=final_q_st, q_en=final_q_en,
                        r_name=ref_name, r_len=ref_len, r_st=final_r_st, r_en=final_r_en,
                        strand='+' if chain_strand else '-',
                        score=score,
                        edit_distance=edit_distance,
                        cigar=cigar
                    )
                    refined_segments.append(segment)

        return refined_segments
        
    def _calculate_edit_distance_from_cigar(self, cigar: str) -> int:
        """Calculate edit distance from a CIGAR string."""
        if not cigar:
            return 0
            
        edit_dist = 0
        num = ""
        for c in cigar:
            if c.isdigit():
                num += c
            else:
                count = int(num)
                num = ""
                # M could be match or mismatch, but we don't have sequence info here
                # so we conservatively count only I, D operations
                if c in 'ID':
                    edit_dist += count
        
        return edit_dist
