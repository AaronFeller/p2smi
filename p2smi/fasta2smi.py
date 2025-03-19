"""
Module to input FASTA peptide files and generate 3D structures.
Original by Fergal; modified by Aaron Feller (2025)
"""

import argparse

import p2smi.utilities.smilesgen as smilesgen


class InvalidConstraintError(Exception):
    pass


def parse_fasta(fasta_file):
    """
    Generator that yields (sequence, constraint) tuples from a FASTA file.
    """
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
    """
    Resolves constraints for a sequence using available functions
    or generates fallback structures.
    """
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
    """
    Processes constraints for all sequences in the FASTA file.
    """
    return (constraint_resolver(seq, constr) for seq, constr in parse_fasta(fasta_file))


def generate_3d_structures(input_fasta, out_file):
    """
    Generates 3D structures from FASTA input and writes to output.
    """
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
    parser = argparse.ArgumentParser(description="Generate peptides from a FASTA file.")
    parser.add_argument(
        "-i", "--input_fasta", required=True, help="FASTA file of peptides."
    )
    parser.add_argument("-o", "--out_file", required=True, help="Output file.")
    args = parser.parse_args()
    generate_3d_structures(args.input_fasta, args.out_file)


if __name__ == "__main__":
    main()
