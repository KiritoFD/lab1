import time
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import os
from collections import defaultdict
import mmap

@dataclass
class RepeatPattern:
    position: int
    length: int
    repeat_count: int
    is_reverse: bool
    original_sequence: str
    query_position: int

class DNARepeatFinder:
    def __init__(self, 
                 similarity_threshold: float = 0.99,
                 min_length: int = 50,
                 max_length: int = 101,
                 max_repeats: int = 1000,
                 prime: int = 31,
                 mod: int = (1 << 31) - 1,
                 base_map: Dict[str, int] = None,
                 query_file: str = "query.txt",
                 reference_file: str = "reference.txt"):
        
        self.similarity_threshold = similarity_threshold
        self.MIN_LENGTH = min_length
        self.MAX_LENGTH = max_length
        self.MAX_REPEATS = max_repeats
        self.prime = prime
        self.mod = mod
        self.query_file = query_file
        self.reference_file = reference_file
        
        self.base_map = base_map or {'A': 3, 'C': 5, 'G': 7, 'T': 11}
        self._precompute_powers(self.MAX_LENGTH)
        
        self.complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        self.rev_comp_cache = {}
        
    def _precompute_powers(self, max_window: int):
        self.powers = [1] * (max_window + 1)
        for i in range(1, max_window + 1):
            self.powers[i] = (self.powers[i-1] * self.prime) % self.mod
            
    def compute_rolling_hashes(self, sequence: str, window_size: int) -> List[int]:
        n = len(sequence)
        if n < window_size:
            return []
        current_hash = 0
        hashes = []
        for i in range(window_size):
            current_hash = (current_hash * self.prime + self.base_map.get(sequence[i], 0)) % self.mod
        hashes.append(current_hash)
        for i in range(window_size, n):
            outgoing_val = self.base_map.get(sequence[i - window_size], 0)
            incoming_val = self.base_map.get(sequence[i], 0)
            current_hash = (current_hash - outgoing_val * self.powers[window_size - 1]) % self.mod
            current_hash = (current_hash * self.prime + incoming_val) % self.mod
            hashes.append(current_hash)
        return hashes

    def build_sequence_hashmap(self, sequence: str, length: int) -> Dict[int, List[int]]:
        hashmap = defaultdict(list)
        hashes = self.compute_rolling_hashes(sequence, length)
        for pos, hash_val in enumerate(hashes):
            hashmap[hash_val].append(pos)
        return hashmap

    def get_reverse_complement(self, seq: str) -> str:
        if seq not in self.rev_comp_cache:
            self.rev_comp_cache[seq] = ''.join(self.complement.get(base, 'N') for base in reversed(seq))
        return self.rev_comp_cache[seq]

    def calculate_similarity(self, seq1: str, seq2: str) -> float:
        return sum(a == b for a, b in zip(seq1, seq2)) / len(seq1)

    def find_consecutive_groups(self, positions: List[int], length: int) -> List[int]:
        if len(positions) < 2:
            return []
        groups = []
        current_count = 1
        for i in range(1, len(positions)):
            if positions[i] == positions[i-1] + length:
                current_count += 1
            else:
                if current_count >= 2:
                    groups.append(current_count)
                current_count = 1
        if current_count >= 2:
            groups.append(current_count)
        return groups

    def filter_nested_repeats(self, repeats: List[RepeatPattern]) -> List[RepeatPattern]:
        to_remove = [False] * len(repeats)
        for i in range(len(repeats)):
            if not to_remove[i]:
                for j in range(i+1, len(repeats)):
                    if (repeats[i].position == repeats[j].position and 
                        repeats[i].is_reverse == repeats[j].is_reverse and
                        repeats[i].length != repeats[j].length):
                        if repeats[i].length > repeats[j].length:
                            to_remove[j] = True
                        else:
                            to_remove[i] = True
                            break
        return [r for i, r in enumerate(repeats) if not to_remove[i]]

    def find_repeats(self, query: str, reference: str) -> List[RepeatPattern]:
        repeats = []
        special_check_around = len(query) // 2
        env_value = os.getenv('SPECIAL_CHECK_AREA')
        if env_value:
            try:
                special_check_around = int(env_value)
            except ValueError:
                pass
        print(f"Special check area around position: {special_check_around}")

        for length in range(self.MIN_LENGTH, min(self.MAX_LENGTH + 1, len(query) + 1)):
            if length > len(reference):
                continue
            query_hashmap = self.build_sequence_hashmap(query, length)
            reference_hashes = self.compute_rolling_hashes(reference, length)
            rev_ref = self.get_reverse_complement(reference)
            rev_ref_hashes = self.compute_rolling_hashes(rev_ref, length)
            reference_len = len(reference)
            
            for i in range(len(reference_hashes)):
                if abs(i - special_check_around) > 10 and length > self.MIN_LENGTH + 10 and length > 100:
                    break

                matches = []
                h_forward = reference_hashes[i]
                if h_forward in query_hashmap:
                    for q_pos in query_hashmap[h_forward]:
                        ref_seg = reference[i:i+length]
                        query_seg = query[q_pos:q_pos+length]
                        if self.calculate_similarity(ref_seg, query_seg) >= self.similarity_threshold:
                            matches.append(q_pos)
                j = reference_len - length - i
                if j >= 0 and j < len(rev_ref_hashes):
                    h_rev_comp = rev_ref_hashes[j]
                    if h_rev_comp in query_hashmap:
                        for q_pos in query_hashmap[h_rev_comp]:
                            rev_comp_seg = self.get_reverse_complement(reference[i:i+length])
                            query_seg = query[q_pos:q_pos+length]
                            if self.calculate_similarity(query_seg, rev_comp_seg) >= self.similarity_threshold:
                                matches.append(q_pos)

                if matches:
                    groups = self.find_consecutive_groups(sorted(matches), length)
                    for group_size in groups:
                        if len(repeats) >= self.MAX_REPEATS:
                            break
                        segment = reference[i:i+length]
                        is_reverse = False
                        repeats.append(RepeatPattern(
                            position=i,
                            length=length,
                            repeat_count=group_size,
                            is_reverse=is_reverse,
                            original_sequence=segment,
                            query_position=matches[0]
                        ))

        repeats = self.filter_nested_repeats(repeats)
        repeats.sort(key=lambda x: x.length * x.repeat_count, reverse=True)
        return repeats[:self.MAX_REPEATS]

def read_sequence(filename: str) -> str:
    with open(filename, 'r') as f:
        try:
            return mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ).read().decode().upper()
        except:
            sequence = f.read().upper()
            return ''.join(filter(str.isalpha, sequence))

def save_repeats_to_file(repeats: List[RepeatPattern]):
    with open('repeat_results_new.txt', 'w') as f:
        f.write("Reference Position,Length,Repeat Count,Is Reverse Repeat,Original Sequence,Query Position\n")
        for r in repeats:
            f.write(f"{r.position},{r.length},{r.repeat_count},"
                   f"{'Yes' if r.is_reverse else 'No'},{r.original_sequence},{r.query_position}\n")

def main(query_file: str = "query.txt", reference_file: str = "reference.txt"):
    print(f"Reading query sequence: {query_file}")
    query = read_sequence(query_file)
    print(f"Reading reference sequence: {reference_file}")
    reference = read_sequence(reference_file)
    
    start_time = time.time()
    finder = DNARepeatFinder(query_file=query_file, reference_file=reference_file)
    repeats = finder.find_repeats(query, reference)
    elapsed_time = time.time() - start_time
    
    print(f"Found {len(repeats)} repeat fragments, elapsed time: {elapsed_time * 1000:.2f} ms")
    
    for i, repeat in enumerate(repeats, 1):
        print(f"Repeat #{i}: Position {repeat.position}, Length {repeat.length}, "
              f"Repeat Count {repeat.repeat_count}, Is Reverse Repeat {'Yes' if repeat.is_reverse else 'No'}")
    
    save_repeats_to_file(repeats)

if __name__ == "__main__":
    import sys
    main(*sys.argv[1:3] if len(sys.argv) >=3 else ())