import sys
from typing import List, Dict, Any

from .data_types import Seed

class Chainer:
    def __init__(self, params: Dict[str, Any]):
        self.params = params # Store gap penalties, thresholds etc.
        # Example params: max_gap, max_dist_diff, chain_score_factor, gap_open, gap_extend

    def _calculate_gap_penalty(self, seed_i: Seed, seed_j: Seed) -> float:
        """Calculates penalty for the gap between two seeds in a chain."""
        q_dist = abs(seed_i.q_pos - seed_j.q_pos)
        r_dist = abs(seed_i.r_pos - seed_j.r_pos)
        dist_diff = abs(q_dist - r_dist)

        # Simple linear gap penalty model (example)
        # A more complex model could use affine gap penalties based on q_dist/r_dist
        # and penalize dist_diff more heavily.
        penalty = (q_dist + r_dist) * self.params.get('gap_penalty_factor', 0.01)
        penalty += dist_diff * self.params.get('dist_diff_penalty_factor', 0.05)

        # Penalize if strands differ (should not happen if check in main loop)
        if seed_i.strand != seed_j.strand:
            penalty += float('inf') # Should not chain across strands

        return penalty

    def chain_seeds(self, seeds: List[Seed]) -> List[List[Seed]]:
        """Chains seeds using dynamic programming."""
        if not seeds:
            return []

        n = len(seeds)
        # Seeds are assumed sorted by q_pos, then ref_id, then r_pos

        dp_score = [self.params.get('k', 15)] * n  # Initial score = kmer match length (approx)
        predecessor = [-1] * n
        max_score_at_i = list(dp_score) # Track overall max ending at i

        max_gap = self.params.get('max_gap', 10000)
        max_dist_diff = self.params.get('max_dist_diff', 500) # Max difference between q_dist and r_dist

        for i in range(n):
            seed_i = seeds[i]
            current_base_score = self.params.get('k', 15) # Base score for the seed itself
            best_prev_score = 0
            best_pred = -1

            # Iterate potential predecessors j < i
            # Optimization: Only check j within a certain query distance window?
            search_start_j = 0 # TODO: Optimize start index based on max_gap?
            for j in range(search_start_j, i):
                seed_j = seeds[j]

                # --- Validity Checks ---
                if seed_i.ref_id != seed_j.ref_id or seed_i.strand != seed_j.strand:
                    continue # Cannot chain across references or strands

                q_dist = seed_i.q_pos - seed_j.q_pos
                r_dist = seed_i.r_pos - seed_j.r_pos

                if not (0 < q_dist < max_gap and 0 < r_dist < max_gap):
                     continue # Must be forward in both, within gap limit

                dist_diff = abs(q_dist - r_dist)
                if dist_diff > max_dist_diff:
                    continue # Check collinearity

                # --- Calculate Score ---
                gap_penalty = self._calculate_gap_penalty(seed_i, seed_j)
                potential_score = max_score_at_i[j] - gap_penalty # Score continues from max chain ending at j

                if potential_score > best_prev_score:
                    best_prev_score = potential_score
                    best_pred = j

            # Update DP table for seed i
            final_score_i = current_base_score + best_prev_score
            dp_score[i] = final_score_i # Score of this specific link, might not be max score ending at i
            max_score_at_i[i] = max(final_score_i, current_base_score) # Max score ending *at or including* i
            if best_prev_score > 0: # Only assign predecessor if we extended a chain
                predecessor[i] = best_pred
            # else: predecessor remains -1


        # --- Extract Chains ---
        chains = []
        visited = [False] * n
        # Sort potential chain ends by score descending to extract best chains first
        potential_ends = sorted(range(n), key=lambda i: max_score_at_i[i], reverse=True)

        min_chain_score = self.params.get('min_chain_score', 50)

        for i in potential_ends:
            if not visited[i] and max_score_at_i[i] >= min_chain_score:
                chain = []
                curr = i
                while curr != -1 and not visited[curr]:
                    visited[curr] = True
                    chain.append(seeds[curr])
                    curr = predecessor[curr]

                if chain: # Should always be true if score > min_chain_score and base score > 0
                    chains.append(chain[::-1]) # Reverse to get correct order

        return chains
