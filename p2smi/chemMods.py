import argparse
import math
import random
import re

from rdkit import Chem


def is_valid_smiles(sequence):
    return Chem.MolFromSmiles(sequence) is not None


def add_n_methylation(sequence, methylation_residue_fraction):
    pattern = r"C\(=O\)N\[C@"
    positions = [m.start() for m in re.finditer(pattern, sequence)]
    num_to_methylate = math.ceil(len(positions) * methylation_residue_fraction)

    for pos in sorted(
        random.sample(positions, min(num_to_methylate, len(positions))),
        reverse=True,
    ):
        sequence = sequence[: pos + 6] + "(C)" + sequence[pos + 6 :]
    return sequence, num_to_methylate


def add_pegylation(sequence):
    peg = "O" + "".join(["CCO" for _ in range(random.randint(1, 4))]) + "C"
    positions = [m.start() for m in re.finditer(r"CN\)", sequence)]
    if not positions:
        return sequence, None
    pos = random.choice(positions)
    return sequence[: pos + 2] + peg + sequence[pos + 2 :], peg


def parse_input_lines(input_lines):
    for line in input_lines:
        parts = line.strip().split(": ", 1)
        if len(parts) == 2:
            yield tuple(parts)
        else:
            yield ("[Malformed line]", None)


def modify_sequence(sequence, methylate, pegylate, nmeth_residues):
    modifications = []
    if methylate:
        sequence, methyl_count = add_n_methylation(sequence, nmeth_residues)
        modifications.append(f"N-methylation({methyl_count})")
    if pegylate:
        sequence, peg = add_pegylation(sequence)
        # find number of PEG units added
        if peg:
            number_of_peg_units = peg.count("CCO")
            modifications.append(f"PEGylation({number_of_peg_units})")
        else:
            modifications.append("PEGylation('N/A')")
    return sequence, modifications


def process_sequences(input_lines, nmeth_rate, peg_rate, nmeth_residues):
    total = len(input_lines)
    methylate_indices = (
        set(random.sample(range(total), math.ceil(total * nmeth_rate)))
        if nmeth_rate > 0
        else set()
    )
    pegylate_indices = (
        set(random.sample(range(total), math.ceil(total * peg_rate)))
        if peg_rate > 0
        else set()
    )

    for i, (header, seq) in enumerate(parse_input_lines(input_lines)):
        if seq is None:
            yield f"{header} [Skipped malformed line]"
            continue

        sequence, mods = modify_sequence(
            seq, i in methylate_indices, i in pegylate_indices, nmeth_residues
        )

        if not is_valid_smiles(sequence):
            yield f"{header} [Invalid SMILES skipped]"
            continue

        mod_str = f" [{' - '.join(mods)}]" if mods else ""
        yield f"{header}{mod_str}: {sequence}"


def process_file(input_file, output_file, peg_rate, nmeth_rate, nmeth_residues):
    with open(input_file, "r") as infile:
        lines = [line.strip() for line in infile if line.strip()]

    modified_lines = process_sequences(lines, nmeth_rate, peg_rate, nmeth_residues)

    with open(output_file, "w") as outfile:
        outfile.write("\n".join(modified_lines) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=("Modify peptide sequences with" "PEGylation and N-methylation.")
    )
    parser.add_argument("-i", "--input_file", required=True, help="Input file path.")
    parser.add_argument("-o", "--output_file", required=True, help="Output file path.")
    parser.add_argument(
        "--peg_rate",
        type=float,
        default=0.2,
        help="Fraction of sequences to PEGylate (0-1).",
    )
    parser.add_argument(
        "--nmeth_rate",
        type=float,
        default=0.2,
        help="Fraction of sequences to N-methylate (0-1).",
    )
    parser.add_argument(
        "--nmeth_residues",
        type=float,
        default=0.2,
        help="Fraction of bonds in each sequence to N-methylate (0-1).",
    )
    args = parser.parse_args()

    process_file(
        args.input_file,
        args.output_file,
        args.peg_rate,
        args.nmeth_rate,
        args.nmeth_residues,
    )


if __name__ == "__main__":
    main()
