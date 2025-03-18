'''
This script generates random amino acid sequences with specified constraints and properties.
It allows the user to specify the length of the sequences, the percentage of noncanonical amino acids,
the percentage of dextro amino acids, and the number of sequences to generate. The generated sequences
can be printed to the console or saved to a file.
'''

import random
import argparse
from p2smi.generate.aminoacids import all_aminos

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
        seq_id = f"seq_{i+1}"
        constraint = random.choice(constraints) if constraints else None
        if constraint:
            seq_id += f"|{constraint}"

        seq_len = random.randint(min_length, max_length)

        # Calculate counts with rounding
        num_dextro_noncanonical = round(seq_len * dextro_percent * noncanonical_percent)
        num_dextro_canonical = round(seq_len * dextro_percent * (1 - noncanonical_percent))
        num_levro_noncanonical = round(seq_len * (1 - dextro_percent) * noncanonical_percent)
        num_levro_canonical = seq_len - (
            num_dextro_noncanonical + num_dextro_canonical + num_levro_noncanonical
        )  # ensure exact sum

        # Build sequence parts
        sequence = []
        sequence.extend(random.choices(canonical_list, k=num_levro_canonical))
        sequence.extend(random.choices(lower_canon_list, k=num_dextro_canonical))
        sequence.extend(random.choices(upper_noncon_list, k=num_levro_noncanonical))
        sequence.extend(random.choices(lower_noncon_list, k=num_dextro_noncanonical))

        random.shuffle(sequence)
        sequences[seq_id] = ''.join(sequence)

    # Print the sequences to stdout or write to a file
    if outfile is None:
        for id, seq in sequences.items():
            print(f">{id}\n{seq}")
    else:
        # Write the sequences to the output file
        with open(outfile, 'w') as f:
            for id, seq in sequences.items():
                f.write(f">{id}\n{seq}\n")            

def main(args):
    # Generate the sequence
    generate_sequence(
        args.min_seq_len, 
        args.max_seq_len, 
        args.noncanonical_percent, 
        args.lowercase_percent,
        args.num_sequences,
        args.constraints,
        args.outfile
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate random amino acid sequences.")
    parser.add_argument("--min_seq_len", type=int, default=5, help="Minimum sequence length (default: 5)")
    parser.add_argument("--max_seq_len", type=int, default=100, help="Maximum sequence length (default: 100)")
    parser.add_argument("--noncanonical_percent", type=float, default=0.2, help="Percentage of noncanonical amino acids (default: 0.2)")
    parser.add_argument("--lowercase_percent", type=float, default=0.2, help="Percentage of dextro amino acids (default: 0.2)")
    parser.add_argument("--num_sequences", type=int, default=1, help="Number of sequences to generate (default: 1)")
    parser.add_argument("--constraints", type=str, nargs='*', help="Constraints for the sequences (use 'all' for all constraints)", default=[])
    parser.add_argument("--outfile", type=str, help="Output file name (default: stdout)", default=None)
    args = parser.parse_args()

    # Handle 'all' option
    if args.constraints == ['all'] or not args.constraints:
        args.constraints = constraint_list

    main(args)