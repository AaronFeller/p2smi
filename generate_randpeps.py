import random
from PepLibGen.StructGen.aminoacids import all_aminos
from tqdm import tqdm
import argparse


def generate_sequence(seq_len, canonical_percent, lowercase_percent):
    # Define the amino acids
    letter_list = []
    for each in all_aminos:
        letter_list+=(all_aminos[each]['Letter'])
    upper_noncon_list = list(set([x.upper() for x in letter_list]))
    canonical_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    upper_noncon_list = [x for x in upper_noncon_list if x not in canonical_list]
    
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
        if i in lowercase_positions and not upper_sequence[i]=='G':
            final_sequence += upper_sequence[i].lower()
        else:
            final_sequence += upper_sequence[i]

    return final_sequence

def main():

    parser = argparse.ArgumentParser(description="Generate random peptide sequences")
    parser.add_argument("-n", "--num_sequences", type=int, help="Number of sequences to generate")
    parser.add_argument("-c", "--canonical_percent", type=float, default=0.9, help="Percentage of canonical amino acids")
    parser.add_argument("-d", "--d_amino_percent", type=float, default=0.1, help="Percentage of lowercase amino acids")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file name")
    args = parser.parse_args()

    
    num_sequences = args.num_sequences
    canonical_percent = args.canonical_percent
    lowercase_percent = args.d_amino_percent
    output_file = args.output

    with open(output_file, 'w') as file:
        for i in tqdm(range(1, num_sequences + 1), desc="Generating sequences"):

            seq_len = random.randint(10, 100)
            sequence = generate_sequence(seq_len, canonical_percent, lowercase_percent)
            cyclization = random.choice(['|HT', '|SS', '|SCNT', '|SCCT', '|SCSC'])
            file.write(">Seq" + str(i) + cyclization + "\n" + sequence + "\n")

if __name__ == "__main__":
    main()
