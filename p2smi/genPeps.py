"""
Random Peptide Sequence Generator

This script generates random amino acid sequences with customizable properties:
- Sequence length (min and max)
- Fraction of noncanonical amino acids
- Fraction of D-amino acids (lowercase residues)
- Optional structural constraints (SS, HT, SCNT, SCCT, SCSC)

Key features:
- Builds sequences using canonical and noncanonical amino acid sets from p2smi.
- Supports randomized constraint assignment or user-specified constraints.
- Outputs sequences in FASTA format to stdout or to a specified output file.

Example usage:
python generate_random_peptides.py \
    --num 10 \
    --min 8 \
    --max 20 \
    --ncaa 0.3 \
    -d 0.25 \
    --constraints all \
    --outfile random_peptides.fasta
"""

import argparse
import random

from p2smi.utilities.aminoacids import all_aminos


def get_amino_acid_lists():
    # Create four lists:
    # 1) canonical amino acids (uppercase)
    # 2) canonical amino acids (lowercase, except Glycine)
    # 3) noncanonical amino acids (uppercase)
    # 4) noncanonical amino acids (lowercase)
    all_aas = [letter for entry in all_aminos.values() for letter in entry["Letter"]]
    canonical = [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    ]
    lower_canon = [aa.lower() for aa in canonical if aa != "G"]
    upper_noncanonical = [
        aa
        for aa in all_aas
        if aa not in canonical and aa not in lower_canon and aa.isupper()
    ]
    lower_noncanonical = [
        aa
        for aa in all_aas
        if aa not in canonical and aa not in lower_canon and aa.islower()
    ]
    return canonical, lower_canon, upper_noncanonical, lower_noncanonical


CONSTRAINTS = ["SS", "HT", "SCNT", "SCCT", "SCSC"]  # Supported constraint types


def calculate_amino_acid_counts(seq_len, noncanonical_percent, dextro_percent):
    # Calculate counts of each category (canonical/noncanonical, L-/D-form)
    dn_nc = round(seq_len * dextro_percent * noncanonical_percent)
    dn_c = round(seq_len * dextro_percent * (1 - noncanonical_percent))
    ln_nc = round(seq_len * (1 - dextro_percent) * noncanonical_percent)
    ln_c = seq_len - (dn_nc + dn_c + ln_nc)
    return ln_c, dn_c, ln_nc, dn_nc


def build_sequence(seq_len, noncanonical_percent, dextro_percent, amino_lists):
    # Build a random sequence of specified length with defined fractions of
    # canonical/noncanonical and D-/L-form residues.
    canonical, lower_canon, upper_noncon, lower_noncon = amino_lists
    ln_c, dn_c, ln_nc, dn_nc = calculate_amino_acid_counts(
        seq_len, noncanonical_percent, dextro_percent
    )
    sequence_parts = (
        random.choices(canonical, k=ln_c)
        + random.choices(lower_canon, k=dn_c)
        + random.choices(upper_noncon, k=ln_nc)
        + random.choices(lower_noncon, k=dn_nc)
    )
    random.shuffle(sequence_parts)
    return "".join(sequence_parts)


def generate_sequences(
    num_sequences,
    min_length,
    max_length,
    noncanonical_percent,
    dextro_percent,
    constraints,
):
    # Generate a dictionary of random sequences with optional constraints
    amino_lists = get_amino_acid_lists()

    def make_sequence(i):
        seq_id = f"seq_{i + 1}"
        constraint = random.choice(constraints) if constraints else None
        if constraint:
            seq_id += f"|{constraint}"
        seq_len = random.randint(min_length, max_length)
        seq = build_sequence(seq_len, noncanonical_percent, dextro_percent, amino_lists)
        return seq_id, seq

    return dict(make_sequence(i) for i in range(num_sequences))


def output_sequences(sequences, outfile=None):
    # Print or write sequences in FASTA format to a file if specified
    lines = [f">{seq_id}\n{seq}" for seq_id, seq in sequences.items()]
    output = "\n".join(lines)
    if outfile:
        with open(outfile, "w") as f:
            f.write(output + "\n")
    else:
        print(output)


def main():
    # CLI entry point for sequence generation with configurable parameters
    parser = argparse.ArgumentParser(
        description="Generate random amino acid sequences."
    )
    parser.add_argument("-n", "--num", type=int, default=10)
    parser.add_argument("-min", "--min_length", type=int, default=10)
    parser.add_argument("-max", "--max_length", type=int, default=10)
    parser.add_argument("-ncaa", "--noncanonical", type=float, default=0)
    parser.add_argument("-d", "--dextro", type=float, default=0)
    parser.add_argument(
        "-c",
        "--cyclization_constraints",
        type=str,
        default=None,
        help="Cyclization types: 'all', 'none', or comma-separated list like 'HT,SCSC'",
    )
    parser.add_argument("-o", "--outfile", type=str, default=None)
    args = parser.parse_args()

    # if constraints is "all", use all supported constraints
    if args.cyclization_constraints == "all":
        constraints = ["SS", "HT", "SCNT", "SCCT", "SCSC"]
    elif args.cyclization_constraints == None:
        constraints = []
    else:
        constraints = [args.cyclization_constraints]

    sequences = generate_sequences(
        args.num,
        args.min_length,
        args.max_length,
        args.noncanonical,
        args.dextro,
        constraints,
    )

    output_sequences(sequences, args.outfile)


if __name__ == "__main__":
    main()
