"""
Module to allow input of fasta formatted files of peptides and generation
of 
Created on 22 Jul 2011

@author: Fergal
"""

import sys

# 3rd party modules
import StructGen
from Bio import SeqIO
class InvalidConstraintError(Exception):
    pass

def parse_fasta(fasta_f):
    """
    Reads a fasta file and returns the sequence.
    """
    with open(fasta_f, "r") as fasta:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            if "|" in seq_record.id:
                constraint = seq_record.id.split("|")[-1]
            else:
                constraint = seq_record.id
            yield str(seq_record.seq), constraint

def process_constraints(fasta_f):
    """
    Processes the fasta file, and fills out incomplete constraint information
    with a best guess.
    """
    error_string = "%s is not a valid constraint for peptide %s"
    constraint_functions = {
        "SS": StructGen.can_ssbond,
        "HT": StructGen.can_htbond,
        "SCNT": StructGen.can_scntbond,
        "SCCT": StructGen.can_scctbond,
        "SCSC": StructGen.can_scscbond,
    }
    for sequence, constraint in parse_fasta(fasta_f):
        # Check if constraint is valid
        if constraint.upper() in StructGen.what_constraints(sequence):
            yield sequence, constraint
        # If constraint does not specify a full pattern, only a type,
        # try and generate the full pattern.
        elif constraint.upper() in list(constraint_functions.keys()):
            result = constraint_functions[constraint.upper()](sequence)
            if result:
                yield result
            else:
                raise InvalidConstraintError(error_string % (sequence, constraint))
        elif constraint.upper() == "SC":
            # Try all the different side chain constraints, return the
            # first one found
            func_keys = [k for k in list(constraint_functions.keys()) if "SC" in k]
            found = False
            for k in func_keys:
                result = constraint_functions[k](sequence)
                if result:
                    yield result
                    found = True
                    break
            if not found:
                raise InvalidConstraintError(error_string % (sequence, constraint))
        elif not constraint:
            # Generate linear peptide
            yield sequence, ""
        else:
            raise InvalidConstraintError(error_string % (sequence, constraint))

def main(fasta_f, out_f="out.sdf"):
    """
    Writes the sequences (with constraints) in fasta_f to a 3d SDF file, out_f.
    """
    out_gen = (
        StructGen.constrained_peptide_smiles(*pair)
        for pair in process_constraints(fasta_f)
    )
    # pdb.set_trace()
    StructGen.write_library(out_gen, out_f, write="structure", write_to_file=True)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate peptides from a fasta file.")
    parser.add_argument("fasta_f", help="Fasta file of peptides.")
    parser.add_argument("--out_f", help="Output file.", default="out.sdf")
    args = parser.parse_args()
    
    main(args.fasta_f, args.out_f)
