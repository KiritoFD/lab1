def read_repeat_details(filename):
    """Reads repeat details from the given file."""
    repeats = []
    with open(filename, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split(',')
            repeats.append({
                'position': int(parts[0]),
                'length': int(parts[1]),
                'repeat_count': int(parts[2]),
                'is_reverse': parts[3] == '是',
                'original_sequence': parts[4],
                'repeat_instance': parts[5]
            })
    return repeats

def read_sequence(filename):
    """Reads a sequence from the given file."""
    with open(filename, 'r') as f:
        return f.read().strip()

if __name__ == "__main__":
    repeat_details_file = "repeat_details.txt"
    reference_file = "reference.txt"
    query_file = "query.txt"

    repeats = read_repeat_details(repeat_details_file)
    reference_sequence = read_sequence(reference_file)
    query_sequence = read_sequence(query_file)

    print(f"Number of repeats found: {len(repeats)}")
    print(f"Length of reference sequence: {len(reference_sequence)}")
    print(f"Length of query sequence: {len(query_sequence)}")

    # Example: Print the first 5 repeats
    print("\nFirst 5 repeats:")
    for i in range(min(5, len(repeats))):
        print(repeats[i])
