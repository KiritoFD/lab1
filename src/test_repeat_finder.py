import os
import random
import time
import subprocess
import csv
from datetime import datetime

def generate_random_dna(length):
    """Generate a random DNA sequence of specified length"""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))

def get_reverse_complement(sequence):
    """Get the reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))

class RepeatRecord:
    def __init__(self, position, sequence, is_reverse=False, count=1):
        self.position = position
        self.sequence = sequence
        self.length = len(sequence)
        self.is_reverse = is_reverse
        self.count = count
    
    def __str__(self):
        return f"Position: {self.position}, Length: {self.length}, Count: {self.count}, " + \
               f"Reverse: {'Yes' if self.is_reverse else 'No'}, Sequence: {self.sequence}"

def generate_test_sequences(ref_length=10000, query_length=10000, num_repeats=50):
    """
    Generate reference and query sequences with deliberate repeats
    Returns: reference sequence, query sequence, list of inserted repeats
    """
    print(f"Generating test sequences: ref={ref_length}bp, query={query_length}bp")
    
    # Generate base sequences
    reference = generate_random_dna(ref_length)
    query = generate_random_dna(query_length)
    
    # Track inserted repeats
    inserted_repeats = []
    
    # Convert to lists for easier manipulation
    ref_list = list(reference)
    query_list = list(query)
    
    # Insert repeats of varying lengths
    repeat_lengths = [5, 10, 15, 20, 30, 50, 75]
    
    for i in range(num_repeats):
        # Choose a random length for this repeat
        length = random.choice(repeat_lengths)
        
        # Extract a segment from reference to use as a repeat
        if ref_length > length:
            ref_pos = random.randint(0, ref_length - length - 1)
            repeat_seq = reference[ref_pos:ref_pos + length]
            
            # Insert into query at a random position
            query_pos = random.randint(0, query_length - length - 1)
            
            # Decide if we'll use reverse complement
            use_reverse = random.choice([True, False])
            
            if use_reverse:
                insert_seq = get_reverse_complement(repeat_seq)
            else:
                insert_seq = repeat_seq
            
            # Insert the sequence
            for j in range(length):
                query_list[query_pos + j] = insert_seq[j]
            
            # Record the repeat
            repeat_record = RepeatRecord(ref_pos, repeat_seq, use_reverse)
            inserted_repeats.append(repeat_record)
            
            print(f"Inserted repeat: {repeat_record}")
    
    # Convert back to strings
    query = ''.join(query_list)
    
    return reference, query, inserted_repeats

def write_sequences_to_files(reference, query, ref_file="reference.txt", query_file="query.txt"):
    """Write sequences to files in FASTA format"""
    with open(ref_file, 'w') as f:
        f.write(f">Reference sequence length {len(reference)}\n")
        for i in range(0, len(reference), 80):
            f.write(reference[i:i+80] + '\n')
    
    with open(query_file, 'w') as f:
        f.write(f">Query sequence length {len(query)}\n")
        for i in range(0, len(query), 80):
            f.write(query[i:i+80] + '\n')

def write_expected_results(repeats, filename="expected_repeats.csv"):
    """Write the expected repeats to a CSV file"""
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Position", "Length", "Count", "ReverseComplement", "Sequence"])
        
        for repeat in repeats:
            writer.writerow([
                repeat.position,
                repeat.length,
                repeat.count,
                "Yes" if repeat.is_reverse else "No",
                repeat.sequence
            ])

def run_repeat_finder(ref_file="reference.txt", query_file="query.txt"):
    """Run the DNA repeat finder program"""
    start_time = time.time()
    
    # Run the program as a subprocess
    print(f"\nRunning DNA repeat finder on {ref_file} and {query_file}...")
    
    # Ensure file paths exist
    if not os.path.exists(ref_file):
        print(f"Error: Reference file not found: {ref_file}")
        return False
    
    if not os.path.exists(query_file):
        print(f"Error: Query file not found: {query_file}")
        return False
    
    try:
        # Use more direct command execution
        cmd = f"dna_repeat_finder.exe \"{ref_file}\" \"{query_file}\""
        print(f"Executing: {cmd}")
        
        result = subprocess.run(
            cmd,
            capture_output=True, 
            text=True,
            shell=True,
            timeout=300  # 5 minute timeout
        )
        
        print(f"Process exited with code: {result.returncode}")
        
        if result.returncode != 0:
            print(f"Error: DNA repeat finder exited with code {result.returncode}")
            print(f"stderr: {result.stderr}")
            return False
        
        print(result.stdout)
    except subprocess.TimeoutExpired:
        print("Error: DNA repeat finder timed out after 5 minutes")
        return False
    except Exception as e:
        print(f"Exception running DNA repeat finder: {e}")
        return False
    
    end_time = time.time()
    print(f"Total execution time: {end_time - start_time:.2f} seconds")
    
    return True

def compare_results(expected_file="expected_repeats.csv", actual_file="repeat_details.txt"):
    """Compare expected results with actual results"""
    print("\nComparing results...")
    
    # Load expected results
    expected_repeats = []
    with open(expected_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            expected_repeats.append(row)
    
    # Load actual results
    actual_repeats = []
    with open(actual_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            actual_repeats.append(row)
    
    # Simple stats
    print(f"Expected repeats: {len(expected_repeats)}")
    print(f"Found repeats: {len(actual_repeats)}")
    
    # More detailed analysis would require matching repeats
    # This would need to consider position shifts, overlaps, etc.
    # For this test, just report the counts
    
    return len(expected_repeats), len(actual_repeats)

def main():
    random.seed(42)  # For reproducibility
    
    # Create test directory
    test_dir = f"test_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    os.makedirs(test_dir, exist_ok=True)
    
    # Test parameters - starting with smallest test first
    test_sizes = [
        (100, 100, 5),    # Tiny test
        (500, 500, 10),   # Very small test
        (1000, 1000, 10), # Small test
    ]
    
    for ref_size, query_size, num_repeats in test_sizes:
        print(f"\n{'='*60}")
        print(f"TEST: Reference={ref_size}bp, Query={query_size}bp, Repeats={num_repeats}")
        print(f"{'='*60}")
        
        # Generate test files
        ref_file = os.path.join(test_dir, f"reference_{ref_size}.txt")
        query_file = os.path.join(test_dir, f"query_{query_size}.txt")
        expected_file = os.path.join(test_dir, f"expected_{ref_size}_{query_size}.csv")
        
        # Generate sequences with known repeats
        reference, query, repeats = generate_test_sequences(ref_size, query_size, num_repeats)
        
        # Write to files
        write_sequences_to_files(reference, query, ref_file, query_file)
        write_expected_results(repeats, expected_file)
        
        # Run repeat finder
        success = run_repeat_finder(ref_file, query_file)
        
        if not success:
            print(f"Test failed for size {ref_size}x{query_size}")
            break
            
        # Compare results if successful
        if success and os.path.exists("repeat_details.txt"):
            expected_count, actual_count = compare_results(expected_file, "repeat_details.txt")
            
            # Calculate detection rate (simple metric)
            if expected_count > 0:
                detection_rate = min(actual_count / expected_count, 1.0) * 100
                print(f"Detection rate: {detection_rate:.1f}%")
            
            # Copy results to test directory for reference
            result_copy = os.path.join(test_dir, f"results_{ref_size}_{query_size}.txt")
            details_copy = os.path.join(test_dir, f"details_{ref_size}_{query_size}.txt")
            
            try:
                with open("repeat_results.txt", 'r') as src, open(result_copy, 'w') as dst:
                    dst.write(src.read())
                with open("repeat_details.txt", 'r') as src, open(details_copy, 'w') as dst:
                    dst.write(src.read())
            except:
                print("Warning: Could not copy result files")

if __name__ == "__main__":
    main()
