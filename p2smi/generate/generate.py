# generate_sequence.py

import random
import argparse
from generate.aminoacids import all_aminos

# Define the canonical and noncanonical amino acids lists
all_aas = []
for each in all_aminos:
    all_aas += (all_aminos[each]['Letter'])

print()

canonical_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
lower_canon_list = [aa.lower() for aa in canonical_list if aa != 'G']

upper_noncon_list = [aa for aa in all_aas if aa not in canonical_list and aa not in lower_canon_list and aa.isupper()]
lower_noncon_list = [aa for aa in all_aas if aa not in canonical_list and aa not in lower_canon_list and aa.islower()]

constraint_list = ['SS', 'HT', 'SCNT', 'SCCT', 'SCSC']


def generate_sequence(min_length, max_length, noncanonical_percent, dextro_percent, num_sequences, constraints, outfile):
    
    sequences = {}

    for i in range(num_sequences):
        # Generate a unique identifier for the sequence        
        id = f"seq_{i+1}"
        
        if constraints != None:
            # If constraints are provided, select one randomly
            constraint = random.choice(args.constraints)
            id += f"|{constraint}"
        
        # Generate a random length for the sequence
        seq_len = random.randint(min_length, max_length)
        
        # Create a list to hold the amino acids
        sequence = []
        
        # Check if the sum of noncanonical and dextro percentages exceeds 1
        if noncanonical_percent + dextro_percent > 1:
            raise ValueError("The sum of noncanonical and dextro percentages cannot exceed 1.")

        # select the amino acids for the sequence
        num_levro_canonical = int(seq_len * (1 - dextro_percent - noncanonical_percent))
        num_dextro_canonical = int(seq_len * (1 - dextro_percent))
        num_levro_noncanonical = int(seq_len * (1 - noncanonical_percent - dextro_percent))
        num_dextro_noncanonical = int(seq_len * noncanonical_percent)
        
        # ensure the sum of all parts equals the total length
        while num_levro_canonical + num_dextro_canonical + num_levro_noncanonical + num_dextro_noncanonical < seq_len:
            num_levro_canonical += 1
        while num_levro_canonical + num_dextro_canonical + num_levro_noncanonical + num_dextro_noncanonical > seq_len:
            num_levro_canonical -= 1
        
        
        # Add the amino acids to the sequence
        if num_levro_canonical > 0:
            sequence.extend(random.choices(canonical_list, k=num_levro_canonical)) 
        if num_dextro_canonical > 0:
            sequence.extend(random.choices(lower_canon_list, k=num_dextro_canonical))
        if num_levro_noncanonical > 0:
            sequence.extend(random.choices(upper_noncon_list, k=num_levro_noncanonical))
        if num_dextro_noncanonical > 0:
            sequence.extend(random.choices(lower_noncon_list, k=num_dextro_noncanonical))
        
        # Shuffle the sequence
        random.shuffle(sequence)
        
        # Join the list into a string
        sequence = ''.join(sequence)
        
        # Store the sequence in the dictionary
        sequences[id] = sequence

    if outfile is None:
        # Print the sequences to stdout
        for id, seq in sequences.items():
            print(f">{id}\n{seq}")
    else:
        # Write the sequences to the output file
        with open(outfile, 'w') as f:
            for id, seq in sequences.items():
                f.write(f">{id}\n{seq}\n")            

def main(args):
    # Generate the sequence
    generate_sequence(args.min_seq_len, 
                        args.max_seq_len, 
                        args.noncanonical_percent, 
                        args.lowercase_percent,
                        args.num_sequences,
                        args.constraints,
                        args.outfile)


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Generate a random amino acid sequence.")
    parser.add_argument("--min_seq_len", type=int, default=5, help="Length of the sequence (default: 5)")
    parser.add_argument("--max_seq_len", type=int, default=100, help="Length of the sequence (default: 100)")
    parser.add_argument("--noncanonical_percent", type=float, default=0.2, help="Percentage of noncanonical amino acids (default: 0.2)")
    parser.add_argument("--lowercase_percent", type=float, default=0.2, help="Percentage of dextro amino acids (default: 0.2)")
    parser.add_argument("--num_sequences", type=int, default=1, help="Number of sequences to generate (default: 1)")
    parser.add_argument("--constraints", type=str, nargs='*', help="Constraints for the sequences (default: None)", choices=constraint_list)
    parser.add_argument("--outfile", type=str, help="Output file name (default: to stdout)", default=None)
    args = parser.parse_args()

    main(args)