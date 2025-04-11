from typing import List, Dict

from .utils import get_minimizers
from .data_types import Seed, ReferenceLocation

class Seeder:
    def __init__(self, ref_index: Dict[int, List[ReferenceLocation]], ref_info: Dict[int, tuple]):
        self.ref_index = ref_index
        self.ref_info = ref_info

    def find_seeds(self, query_seq: str, k: int, w: int) -> List[Seed]:
        """Finds seeds (minimizer matches) between query and reference index."""
        seeds = []
        query_minimizers = get_minimizers(query_seq, k, w)

        for q_mini_hash, q_pos, q_is_forward in query_minimizers:
            if q_mini_hash in self.ref_index:
                for ref_loc in self.ref_index[q_mini_hash]:
                    # Determine relative strand:
                    # Match if query strand matches ref strand (both True or both False)
                    same_strand = (q_is_forward == ref_loc.strand)
                    
                    # Get reference name for this seed
                    ref_name = self.ref_info[ref_loc.ref_id][0] if ref_loc.ref_id in self.ref_info else "unknown"
                    
                    seed = Seed(
                        q_pos=q_pos,
                        ref_id=ref_loc.ref_id,
                        r_pos=ref_loc.pos,
                        strand=same_strand,  # True if same relative strand, False if opposite
                        ref_name=ref_name
                    )
                    seeds.append(seed)

        # Sort by query position, then ref position (using Seed.__lt__)
        seeds.sort() 
        return seeds
