import random
import math
import re
import argparse


def add_n_methylation(sequence, methylation_percent=0.2):
    """
    Adds N-methylation to a percentage of peptidic bonds in a sequence.

    Args:
        sequence (str): The input peptide sequence in SMILES format.
        methylation_percent (float): Percentage of peptidic bonds to methylate (default is 20%).

    Returns:
        str: The methylated peptide sequence.
        int: Number of N-methylations applied.
    """
    num_peptidic_bonds = sequence.count('C(=O)N[C@')
    num_methylate = math.ceil(num_peptidic_bonds * methylation_percent)
    peptidic_bond_positions = [m.start() for m in re.finditer(r'C\(=O\)N\[C@', sequence)]

    methylated_seq = sequence
    if num_methylate > 0 and peptidic_bond_positions:
        methylated_positions = random.sample(peptidic_bond_positions, num_methylate)
        for position in methylated_positions:
            methylated_seq = (
                methylated_seq[:position + 6] + '(C)' + methylated_seq[position + 6:]
            )

    return methylated_seq, num_methylate


def add_pegylation(sequence):
    """
    Adds PEGylation to one free amine group in a sequence.

    Args:
        sequence (str): The input peptide sequence in SMILES format.

    Returns:
        str: The PEGylated peptide sequence.
        str: The PEG fragment used.
    """
    PEG_fragment = 'O' + random.randint(1, 4) * 'CCO' + 'C'
    free_amine_positions = [m.start() for m in re.finditer(r'CN\)', sequence)]

    if not free_amine_positions:
        raise ValueError("No free amine groups (CN) found for PEGylation.")

    pegylated_position = random.choice(free_amine_positions)
    pegylated_seq = (
        sequence[:pegylated_position + 2] + PEG_fragment + sequence[pegylated_position + 2:]
    )

    return pegylated_seq, PEG_fragment


def process_file(input_file, output_file, methylation_percent, apply_methylation, apply_pegylation):
    """
    Processes a file containing peptide sequences in SMILES format, applies modifications,
    and writes the modified sequences to a new file.

    Args:
        input_file (str): Path to the input file.
        output_file (str): Path to the output file.
        methylation_percent (float): Percentage for N-methylation.
        apply_methylation (bool): Whether to apply N-methylation.
        apply_pegylation (bool): Whether to apply PEGylation.
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            # Skip empty lines
            if not line.strip():
                continue

            # Split the header and sequence
            header, sequence = line.strip().split(": ", 1)

            modifications = []

            # Apply N-methylation
            if apply_methylation:
                sequence, num_methylations = add_n_methylation(sequence, methylation_percent)
                modifications.append(f"N-methylation({num_methylations})")

            # Apply PEGylation
            if apply_pegylation:
                try:
                    sequence, peg_fragment = add_pegylation(sequence)
                    modifications.append(f"PEGylation({peg_fragment})")
                except ValueError as e:
                    modifications.append("PEGylation(N/A)")

            # Update the header with modifications
            modified_header = f"{header} [{'-'.join(modifications)}]"

            # Write the modified sequence to the output file
            outfile.write(f"{modified_header}: {sequence}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modify peptide sequences in a SMILES file.")
    parser.add_argument("-i", "--input_file", required=True, help="Path to the input file.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output file.")
    parser.add_argument("--n_methylation", action="store_true", help="Apply N-methylation to peptidic bonds.")
    parser.add_argument("--pegylation", action="store_true", help="Apply PEGylation to a free amine group.")
    parser.add_argument("--percent", type=float, default=0.2, help="Percentage for N-methylation (default: 20%).")
    args = parser.parse_args()
    
    
    if not args.n_methylation and not args.pegylation:
        raise ValueError("At least one modification type must be specified (--n_methylation or --pegylation).")

    process_file(
        args.input_file,
        args.output_file,
        args.percent,
        args.n_methylation,
        args.pegylation
    )