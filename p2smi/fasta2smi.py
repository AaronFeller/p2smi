"""
Module to input FASTA peptide files and generate 3D structures.

Original by Fergal; modified by Aaron Feller (2025).
Reads peptide sequences from a FASTA file, resolves structural constraints,
converts them to SMILES strings, and writes outputs.

Uses p2smi.utilities.smilesgen.
"""

import argparse
import p2smi.utilities.smilesgen as smilesgen


class InvalidConstraintError(Exception):
    # Custom exception for invalid constraints
    pass


def parse_fasta(fasta_file):
    # Parse a FASTA file and yield (sequence, constraint) tuples.
    # Constraint is taken from the header line after a '|' if present.
    with open(fasta_file, "r") as fasta:
        sequence, constraint = "", ""
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    yield sequence, constraint
                sequence = ""
                constraint = line.split("|")[-1] if "|" in line else ""
            else:
                sequence += line
        if sequence:
            yield sequence, constraint


def constraint_resolver(sequence, constraint):
    # Resolve constraints by checking known patterns or fallback attempts.
    # Return (sequence, constraint) or fallback to linear if none apply.
    constraint_functions = {
        "SS": smilesgen.can_ssbond,
        "HT": smilesgen.can_htbond,
        "SCNT": smilesgen.can_scntbond,
        "SCCT": smilesgen.can_scctbond,
        "SCSC": smilesgen.can_scscbond,
    }
    valid_constraints = smilesgen.what_constraints(sequence)

    if constraint.upper() in valid_constraints:
        return (sequence, constraint)
    elif constraint.upper() in constraint_functions:
        result = constraint_functions[constraint.upper()](sequence)
        return result or (sequence, "")
    elif constraint.upper() == "SC":
        # If "SC" is provided, try each SC-related constraint function
        for func in [
            constraint_functions[k] for k in constraint_functions if "SC" in k
        ]:
            result = func(sequence)
            if result:
                return result
        raise InvalidConstraintError(f"{sequence} has invalid constraint {constraint}")
    elif constraint in (None, ""):
        return (sequence, "")
    else:
        raise InvalidConstraintError(f"{sequence} has invalid constraint {constraint}")


def process_constraints(fasta_file):
    # Process all sequences from the FASTA file through constraint resolution
    return (constraint_resolver(seq, constr) for seq, constr in parse_fasta(fasta_file))


def generate_smiles_strings(input_fasta, out_file):
    # Generate SMILES and write peptide structures from FASTA input to output file
    resolved_sequences = process_constraints(input_fasta)
    smilesgen.write_library(
        (
            smilesgen.constrained_peptide_smiles(seq, constr)
            for seq, constr in resolved_sequences
        ),
        out_file,
        write="text",
        write_to_file=True,
    )


def main():
    # CLI entry point: takes FASTA file input, output file path, and generates structures
    parser = argparse.ArgumentParser(
        description="Generate peptides from a standard FASTA file."
    )
    parser.add_argument(
        "-i", "--input_fasta", required=True, help="FASTA file of peptides."
    )
    parser.add_argument("-o", "--out_file", required=True, help="Output file.")
    args = parser.parse_args()

    generate_smiles_strings(args.input_fasta, args.out_file)


if __name__ == "__main__":
    main()
