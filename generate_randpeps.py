import random

def generate_sequence(seq_len, canonical_percent, lowercase_percent):
    # Define the amino acids
    canonical_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    upper_noncon_list = ['B', 'J', 'O', 'U', 'X', 'Z']

    # Calculate the percentage of noncanonical amino acids
    noncanonical_percent = 1 - canonical_percent

    # Define the number of canonical and noncanonical amino acids
    num_canonical = int(seq_len * canonical_percent)
    num_noncanonical = seq_len - num_canonical

    # Generate a list of canonical and noncanonical amino acids
    canonical_amino_acids = random.choices(canonical_list, k=num_canonical)
    noncanonical_amino_acids = random.choices(upper_noncon_list, k=num_noncanonical)

    # Combine the canonical and noncanonical amino acids
    sequence = canonical_amino_acids + noncanonical_amino_acids

    # Shuffle the sequence
    random.shuffle(sequence)

    # Join the amino acids to form the final sequence
    upper_sequence = ''.join(sequence)

    # Make a percentage of the sequence lowercase
    lowercase_positions = random.sample(range(seq_len), int(seq_len * lowercase_percent))
    final_sequence = ''
    for i in range(seq_len):
        if i in lowercase_positions:
            final_sequence += upper_sequence[i].lower()
        else:
            final_sequence += upper_sequence[i]

    return final_sequence

def main(args):
    num_sequences = args.num_sequences
    canonical_percent = args.canonical_percent
    lowercase_percent = args.d_amino_percent
    output_file = args.output

    with open(output_file, 'w') as file:
        for i in range(1, num_sequences + 1):
            seq_len = random.randint(1, 100)
            sequence = generate_sequence(seq_len, canonical_percent, lowercase_percent)
            file.write(f">Seq{i}\n{sequence}\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate random peptide sequences")
    parser.add_argument("num_sequences", type=int, help="Number of sequences to generate")
    parser.add_argument("--canonical_percent", type=float, default=0.9, help="Percentage of canonical amino acids")
    parser.add_argument("--d_amino_percent", type=float, default=0.1, help="Percentage of lowercase amino acids")
    parser.add_argument("--output", type=str, default='output.fasta', help="Output file name")
    args = parser.parse_args()
    
    main(args)
