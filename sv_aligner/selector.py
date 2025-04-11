import bisect
from typing import List

from .data_types import AlignmentSegment

class FinalSelector:

    def select_alignments(self, segments: List[AlignmentSegment]) -> List[AlignmentSegment]:
        """Selects the best non-overlapping set of alignments based on score."""
        if not segments:
            return []

        # Sort segments by query end position
        segments.sort(key=lambda s: s.q_en)
        n = len(segments)

        # Find the latest non-overlapping predecessor for each segment
        # p[i] = index j < i such that segments[j].q_en <= segments[i].q_st (max j)
        predecessor = [-1] * n
        # Create list of end points for efficient search
        end_points = [s.q_en for s in segments]

        for i in range(n):
            # Find the rightmost index j such that segments[j].q_en <= segments[i].q_st
            # bisect_right finds insertion point for segments[i].q_st in end_points
            # The element to the left (index - 1) is the one we want.
            insertion_point = bisect.bisect_right(end_points, segments[i].q_st, hi=i) # Search only up to i
            if insertion_point > 0:
                predecessor[i] = insertion_point - 1
            # else: predecessor remains -1

        # Dynamic programming to find max score
        max_score = [0.0] * n
        # Store which decision was taken: 0=exclude segment i, 1=include segment i
        decision = [0] * n

        for i in range(n):
            score_i = segments[i].score # Use the alignment score
            prev_max_score = max_score[predecessor[i]] if predecessor[i] != -1 else 0
            score_including_i = score_i + prev_max_score

            score_excluding_i = max_score[i-1] if i > 0 else 0

            if score_including_i > score_excluding_i:
                max_score[i] = score_including_i
                decision[i] = 1 # Include segment i
            else:
                max_score[i] = score_excluding_i
                decision[i] = 0 # Exclude segment i

        # Backtrack to reconstruct the optimal set
        selected_segments = []
        i = n - 1
        while i >= 0:
            if decision[i] == 1:
                selected_segments.append(segments[i])
                i = predecessor[i] # Jump to the previous segment in the optimal chain
            else:
                i -= 1 # Move to the next segment without including this one

        return selected_segments[::-1] # Reverse to get original sort order by position
