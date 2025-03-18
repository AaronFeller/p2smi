"""
Module to allow input of fasta formatted files of peptides and generation of 3D structures
Created on 22 Jul 2011

@author: Fergal

Modified by: Aaron Feller (2025)
"""

import smilesgen as smilesgen
import argparse

class InvalidConstraintError(Exception):
    pass

def parse_fasta(fasta_f):
    """
    Reads a fasta file and returns the sequence along with its identifier.
    """
    with open(fasta_f, "r") as fasta:
        sequence = ""
        constraint = ""
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                if sequence:  # Yield previous sequence and identifier
                    yield sequence, constraint
                sequence = ""
                # if line doesn't contain a constraint, leave it as empty string
                if line.startswith(">") and "|" in line:
                    constraint = line.split("|")[-1]
            else:
                sequence += line
        if sequence:  # Yield the last sequence and identifier
            yield sequence, constraint

def process_constraints(fasta_f):
    """
    Processes the fasta file, and fills out incomplete constraint information
    with a best guess.
    """
    error_string = "%s is not a valid constraint for peptide %s"
    constraint_functions = {
        "SS": smilesgen.can_ssbond,
        "HT": smilesgen.can_htbond,
        "SCNT": smilesgen.can_scntbond,
        "SCCT": smilesgen.can_scctbond,
        "SCSC": smilesgen.can_scscbond,
    }
    for sequence, constraint in parse_fasta(fasta_f):
        # Check if constraint is valid
        if constraint.upper() in smilesgen.what_constraints(sequence):
            yield sequence, constraint
        # If constraint does not specify a full pattern, only a type,
        # try and generate the full pattern.
        elif constraint.upper() in list(constraint_functions.keys()):
            result = constraint_functions[constraint.upper()](sequence)
            if result:
                yield result
            else:
                yield sequence, ""
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
        elif constraint == None:
            yield sequence, ""
        else:
            raise InvalidConstraintError(error_string % (sequence, constraint))

def main(args):
    """
    Writes the sequences (with constraints) in fasta_f to a 3d SDF file, out_f.
    """
    fasta = args.input_fasta
    outfile = args.out_file
    
    out_gen = (
        smilesgen.constrained_peptide_smiles(*pair) for pair in process_constraints(fasta)
    )
    
    smilesgen.write_library(out_gen, outfile, write="text", write_to_file=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate peptides from a fasta file.")
    parser.add_argument("-i", "--input_fasta", help="Fasta file of peptides.", required=True)
    parser.add_argument("-o", "--out_file", help="Output file.", required=True)
    
    args = parser.parse_args()
    
    main(args)

