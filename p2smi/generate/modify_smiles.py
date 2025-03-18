import random
import math
import re
import argparse
from rdkit import Chem


def is_valid_smiles(sequence):
    """Checks if the given sequence is a valid SMILES string."""
    mol = Chem.MolFromSmiles(sequence)
    return mol is not None

def add_n_methylation(sequence, methylation_percent=0.2):
    """Adds N-methylation to a percentage of peptidic bonds in a sequence."""
    num_peptidic_bonds = len(re.findall(r'C\(=O\)N\[C@', sequence))  
    num_methylate = math.ceil(num_peptidic_bonds * methylation_percent)
    peptidic_bond_positions = [m.start() for m in re.finditer(r'C\(=O\)N\[C@', sequence)]
    
    if num_methylate > 0 and peptidic_bond_positions:
        methylated_positions = random.sample(peptidic_bond_positions, min(num_methylate, len(peptidic_bond_positions)))
        for position in sorted(methylated_positions, reverse=True):
            sequence = sequence[:position + 6] + '(C)' + sequence[position + 6:]  # Insert methylation

    return sequence, num_methylate


def add_pegylation(sequence):
    """Adds PEGylation to one free amine group in a sequence."""
    PEG_fragment = 'O' + ''.join(['CCO' for _ in range(random.randint(1, 4))]) + 'C'
    free_amine_positions = [m.start() for m in re.finditer(r'CN\)', sequence)]

    if not free_amine_positions:
        return sequence, None  # No PEGylation possible

    pegylated_position = random.choice(free_amine_positions)
    pegylated_seq = sequence[:pegylated_position + 2] + PEG_fragment + sequence[pegylated_position + 2:]

    return pegylated_seq, PEG_fragment


def process_file(input_file, output_file, methylation_percent, apply_methylation, apply_pegylation, methylate_fraction, pegylate_fraction):
    """Processes peptide sequences, applies modifications randomly based on set fractions, and writes output."""

    print(f"Processing file: {input_file}")
    
    with open(input_file, "r") as infile:
        lines = [line.strip() for line in infile if line.strip()]

    total_sequences = len(lines)
    num_methylate_sequences = math.ceil(total_sequences * methylate_fraction)
    num_pegylate_sequences = math.ceil(total_sequences * pegylate_fraction)

    # Randomly select sequences to modify
    methylate_indices = set(random.sample(range(total_sequences), num_methylate_sequences)) if apply_methylation else set()
    pegylate_indices = set(random.sample(range(total_sequences), num_pegylate_sequences)) if apply_pegylation else set()

    with open(output_file, "w") as outfile:
        for i, line in enumerate(lines):
            parts = line.split(": ", 1)
            if len(parts) != 2:
                print(f"[Warning] Skipping malformed line: {line}")
                continue

            header, sequence = parts
            modifications = []

            # Apply N-methylation if this sequence is selected
            if apply_methylation and i in methylate_indices:
                sequence, num_methylations = add_n_methylation(sequence, methylation_percent)
                modifications.append(f"N-methylation({num_methylations})")

            # Apply PEGylation if this sequence is selected
            if apply_pegylation and i in pegylate_indices:
                sequence, peg_fragment = add_pegylation(sequence)
                if peg_fragment:
                    modifications.append(f"PEGylation({peg_fragment})")
                else:
                    modifications.append("PEGylation(N/A)")

            
            # check sequence with RDkit
            if not is_valid_smiles(sequence):
                print(f"[Warning] Invalid SMILES sequence: {sequence}")
                continue
            
            # Update header with modifications
            modified_header = f"{header} [{'-'.join(modifications)}]" if modifications else header

            # Write to output file
            outfile.write(f"{modified_header}: {sequence}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modify peptide sequences in a SMILES file.")
    parser.add_argument("-i", "--input_file", required=True, help="Path to the input file.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output file.")
    parser.add_argument("--n_methylation", action="store_true", help="Apply N-methylation to peptidic bonds.")
    parser.add_argument("--pegylation", action="store_true", help="Apply PEGylation to a free amine group.")
    parser.add_argument("--percent", type=float, default=0.2, help="Percentage for N-methylation (default: 20%).")
    parser.add_argument("--methylate_fraction", type=float, default=0.2, help="Fraction of sequences to methylate (default: 100%).")
    parser.add_argument("--pegylate_fraction", type=float, default=0.2, help="Fraction of sequences to PEGylate (default: 100%).")

    args = parser.parse_args()

    if not args.n_methylation and not args.pegylation:
        raise ValueError("At least one modification type must be specified (--n_methylation or --pegylation).")

    process_file(
        args.input_file,
        args.output_file,
        args.percent,
        args.n_methylation,
        args.pegylation,
        args.methylate_fraction,
        args.pegylate_fraction
    )