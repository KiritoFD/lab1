import random
import time
import argparse

def generate_dna_sequence(length):
    """Generate a random DNA sequence of specified length"""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))

class RepeatInfo:
    def __init__(self, position, repeat_seq, is_reverse_complement=False, count=1):
        self.position = position
        self.sequence = repeat_seq
        self.length = len(repeat_seq)
        self.is_reverse_complement = is_reverse_complement
        self.count = count

def get_reverse_complement(sequence):
    """Get the reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))

def insert_repeats(sequence, min_repeat_len=5, max_repeat_len=120):
    """Insert various types of repeats into the sequence and return tracking information"""
    sequence_list = list(sequence)
    seq_length = len(sequence)
    
    # Number of repeats scales with sequence length
    num_repeats = seq_length // 10000  # One repeat per 10k bases
    
    print(f"Inserting {num_repeats} repeat patterns...")
    
    # Track all inserted repeats
    repeats = []
    
    # Short repeats (5-15 bp)
    for _ in range(num_repeats):
        repeat_len = random.randint(5, 15)
        repeat = ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(repeat_len))
        
        # Insert at 2-4 locations
        num_locations = random.randint(2, 4)
        positions = []
        
        for _ in range(num_locations):
            pos = random.randint(0, seq_length - repeat_len - 1)
            positions.append(pos)
            for j in range(repeat_len):
                sequence_list[pos + j] = repeat[j]
        
        # Add to tracking
        repeats.append(RepeatInfo(min(positions), repeat, False, num_locations))
    
    # Medium repeats (20-50 bp)
    for _ in range(num_repeats // 2):
        repeat_len = random.randint(20, 50)
        repeat = ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(repeat_len))
        
        # Insert at 1-3 locations
        num_locations = random.randint(1, 3)
        positions = []
        
        for _ in range(num_locations):
            pos = random.randint(0, seq_length - repeat_len - 1)
            positions.append(pos)
            for j in range(repeat_len):
                sequence_list[pos + j] = repeat[j]
        
        # Add to tracking
        repeats.append(RepeatInfo(min(positions), repeat, False, num_locations))
    
    # Long repeats (80-120 bp)
    for _ in range(num_repeats // 5):
        repeat_len = random.randint(80, 120)
        repeat = ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(repeat_len))
        
        # Insert at 1-2 locations
        num_locations = random.randint(1, 2)
        positions = []
        
        for _ in range(num_locations):
            pos = random.randint(0, seq_length - repeat_len - 1)
            positions.append(pos)
            for j in range(repeat_len):
                sequence_list[pos + j] = repeat[j]
        
        # Add to tracking
        repeats.append(RepeatInfo(min(positions), repeat, False, num_locations))
    
    # Insert some reverse complement repeats
    for _ in range(num_repeats // 2):
        repeat_len = random.randint(10, 50)
        repeat = ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(repeat_len))
        rev_comp = get_reverse_complement(repeat)
        
        # Insert original and reverse complement
        pos1 = random.randint(0, seq_length // 2 - repeat_len - 1)
        pos2 = (seq_length // 2) + random.randint(0, seq_length // 2 - repeat_len - 1)
        
        for j in range(repeat_len):
            sequence_list[pos1 + j] = repeat[j]
            sequence_list[pos2 + j] = rev_comp[j]
        
        # Add to tracking
        repeats.append(RepeatInfo(pos1, repeat, True, 1))
    
    return ''.join(sequence_list), repeats

def write_sequence_to_file(filename, sequence, header):
    """Write a DNA sequence to a file in FASTA format"""
    with open(filename, 'w') as f:
        f.write(f">{header}\n")
        
        # Write in chunks of 80 characters per line
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + '\n')

def write_answer_key(filename, repeats, ref_repeats=None):
    """Write the answer key with information about inserted repeats"""
    with open(filename, 'w') as f:
        f.write("Position,Length,RepeatCount,ReverseComplement,Sequence\n")
        
        for repeat in repeats:
            f.write(f"{repeat.position},{repeat.length},{repeat.count}," + 
                    f"{'Yes' if repeat.is_reverse_complement else 'No'}," +
                    f"{repeat.sequence}\n")
        
        if ref_repeats:
            f.write("\n# Reference repeats\n")
            for repeat in ref_repeats:
                f.write(f"{repeat.position},{repeat.length},{repeat.count}," + 
                        f"{'Yes' if repeat.is_reverse_complement else 'No'}," +
                        f"{repeat.sequence}\n")

def main():
    parser = argparse.ArgumentParser(description='Generate DNA sequences for testing')
    parser.add_argument('-ref_length', type=int, default=100000, help='Length of reference sequence')
    parser.add_argument('-query_length', type=int, default=100000, help='Length of query sequence')
    parser.add_argument('-ref', default="reference.txt", help='Reference sequence output filename')
    parser.add_argument('-query', default="query.txt", help='Query sequence output filename')
    parser.add_argument('-answers', default="repeat_answers.txt", help='Answer key filename')
    
    args = parser.parse_args()
    
    ref_length = args.ref_length
    query_length = args.query_length
    ref_filename = args.ref
    query_filename = args.query
    answers_filename = args.answers
    
    print(f"Generating DNA sequences: reference={ref_length}bp, query={query_length}bp")
    start_time = time.time()
    
    # Generate reference sequence
    print("Generating reference sequence...")
    reference = generate_dna_sequence(ref_length)
    
    # Generate query sequence
    print("Generating query sequence...")
    query = generate_dna_sequence(query_length)
    
    # Insert repeats in sequences
    print("Inserting repeats in reference sequence...")
    reference, ref_repeats = insert_repeats(reference)
    
    print("Inserting repeats in query sequence...")
    query, query_repeats = insert_repeats(query)
    
    # Add some shared repeats between reference and query
    print("Adding shared repeats between sequences...")
    shared_repeat_count = min(len(ref_repeats), len(query_repeats)) // 3
    
    for i in range(shared_repeat_count):
        # Pick a random repeat from reference
        ref_repeat = ref_repeats[random.randint(0, len(ref_repeats)-1)]
        repeat_seq = ref_repeat.sequence
        
        # Insert it into query at a random position
        pos = random.randint(0, query_length - ref_repeat.length - 1)
        query_list = list(query)
        for j in range(ref_repeat.length):
            query_list[pos + j] = repeat_seq[j]
        query = ''.join(query_list)
        
        # Add to tracking
        query_repeats.append(RepeatInfo(pos, repeat_seq, ref_repeat.is_reverse_complement))
    
    # Add some mutations to differentiate sequences
    if ref_length == query_length:
        print("Adding mutations between sequences...")
        num_differences = ref_length // 20  # About 5% difference
        query_list = list(query)
        for _ in range(num_differences):
            pos = random.randint(0, query_length - 1)
            current_base = query_list[pos]
            new_base = random.choice([b for b in 'ACGT' if b != current_base])
            query_list[pos] = new_base
        query = ''.join(query_list)
    
    end_time = time.time()
    generation_time = end_time - start_time
    
    # Save sequences to files
    print(f"Writing reference sequence to {ref_filename}...")
    write_sequence_to_file(ref_filename, reference, f"Reference sequence length {ref_length}")
    
    print(f"Writing query sequence to {query_filename}...")
    write_sequence_to_file(query_filename, query, f"Query sequence length {query_length}")
    
    # Save answer key
    print(f"Writing repeat answer key to {answers_filename}...")
    write_answer_key(answers_filename, query_repeats, ref_repeats)
    
    print("\nGenerated sequences successfully:")
    print(f"- Reference sequence saved to: {ref_filename}")
    print(f"- Query sequence saved to: {query_filename}")
    print(f"- Answer key saved to: {answers_filename}")
    print(f"- Reference length: {ref_length} base pairs")
    print(f"- Query length: {query_length} base pairs")
    print(f"- Generation time: {generation_time:.2f} seconds")
    print(f"- Total reference repeats: {len(ref_repeats)}")
    print(f"- Total query repeats: {len(query_repeats)}")

if __name__ == "__main__":
    random.seed(time.time())
    main()
