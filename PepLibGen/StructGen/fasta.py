"""
Module to allow input of fasta formatted files of peptides and generation of 3D structures
Created on 22 Jul 2011

@author: Fergal
"""

import sys
import StructGen
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
                constraint = line.split("|")[-1] if "|" in line else line[1:]
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
        else:
            raise InvalidConstraintError(error_string % (sequence, constraint))

# def process_constraints(fasta_f):
#     """
#     Processes the fasta file and fills out incomplete constraint information
#     with a best guess.
#     """
#     error_string = "{} is not a valid constraint for peptide {}"
#     constraint_functions = {
#         "SS": StructGen.can_ssbond,
#         "HT": StructGen.can_htbond,
#         "SCNT": StructGen.can_scntbond,
#         "SCCT": StructGen.can_scctbond,
#         "SCSC": StructGen.can_scscbond,
#     }
    
#     for sequence, constraint in parse_fasta(fasta_f):
#         upper_constraint = constraint.upper()
#         valid_constraints = StructGen.what_constraints(sequence)
#         if upper_constraint in valid_constraints:
#             yield sequence, constraint
#         elif upper_constraint in constraint_functions:
#             if constraint_functions[upper_constraint](sequence):
#                 yield sequence, upper_constraint
#             else:
#                 raise InvalidConstraintError(error_string.format(constraint, sequence))
#         elif upper_constraint == "SC":
#             found = False
#             for k in ["SCNT", "SCCT", "SCSC"]:
#                 if constraint_functions[k](sequence):
#                     yield sequence, k
#                     found = True
#                     break
#             if not found:
#                 raise InvalidConstraintError(error_string.format(constraint, sequence))
#         elif not constraint:
#             yield sequence, "L"  # Assume 'L' for linear if no constraint provided
#         else:
#             raise InvalidConstraintError(error_string.format(constraint, sequence))


def main(args):
    """
    Writes the sequences (with constraints) in fasta_f to a 3d SDF file, out_f.
    """
    fasta = args.input_fasta
    outfile = args.out_file
    
    out_gen = (
        StructGen.constrained_peptide_smiles(*pair) for pair in process_constraints(fasta)
    )
    
    StructGen.write_library(out_gen, outfile, write="text", write_to_file=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate peptides from a fasta file.")
    parser.add_argument("-i", "--input_fasta", help="Fasta file of peptides.", required=True)
    parser.add_argument("-o", "--out_file", help="Output file.", required=True)
    
    args = parser.parse_args()
    
    main(args)


# """
# Module to allow input of fasta formatted files of peptides and generation
# of 
# Created on 22 Jul 2011

# @author: Fergal
# """

# import sys

# # 3rd party modules
# import StructGen
# from Bio import SeqIO
# class InvalidConstraintError(Exception):
#     pass

# def parse_fasta(fasta_f):
#     """
#     Reads a fasta file and returns the sequence along with its identifier.
#     """
#     with open(fasta_f, "r") as fasta:
#         sequence = ""
#         constraint = ""
#         for line in fasta:
#             line = line.strip()
#             if line.startswith(">"):
#                 if sequence:  # Yield previous sequence and identifier
#                     yield sequence, constraint
#                 sequence = ""
#                 if "|" in line:
#                     constraint = line.split("|")[-1]
#                 else:
#                     constraint = line[1:]  # Remove '>' character
#             else:
#                 sequence += line
#         if sequence:  # Yield the last sequence and identifier
#             yield sequence, constraint

# # Example usage:
# # for seq, id in parse_fasta("your_fasta_file.fasta"):
# #     print(f"ID: {id}, Sequence: {seq}")

# def process_constraints(fasta_f):
#     """
#     Processes the fasta file, and fills out incomplete constraint information
#     with a best guess.
#     """
#     error_string = "%s is not a valid constraint for peptide %s"
#     constraint_functions = {
#         "SS": StructGen.can_ssbond,
#         "HT": StructGen.can_htbond,
#         "SCNT": StructGen.can_scntbond,
#         "SCCT": StructGen.can_scctbond,
#         "SCSC": StructGen.can_scscbond,
#     }
#     for sequence, constraint in parse_fasta(fasta_f):
#         # Check if constraint is valid
#         if constraint.upper() in StructGen.what_constraints(sequence, constraint):
#             yield sequence, constraint
#         # If constraint does not specify a full pattern, only a type,
#         # try and generate the full pattern.
#         elif constraint.upper() in list(constraint_functions.keys()):
#             result = constraint_functions[constraint.upper()](sequence)
#             if result:
#                 yield result
#             else:
#                 raise InvalidConstraintError(error_string % (sequence, constraint))
#         elif constraint.upper() == "SC":
#             # Try all the different side chain constraints, return the
#             # first one found
#             func_keys = [k for k in list(constraint_functions.keys()) if "SC" in k]
#             found = False
#             for k in func_keys:
#                 result = constraint_functions[k](sequence)
#                 if result:
#                     yield result
#                     found = True
#                     break
#             if not found:
#                 raise InvalidConstraintError(error_string % (sequence, constraint))
#         elif not constraint:
#             # Generate linear peptide
#             yield sequence, ""
#         else:
#             raise InvalidConstraintError(error_string % (sequence, constraint))

# def main(fasta_f, out_f="out.sdf"):
#     """
#     Writes the sequences (with constraints) in fasta_f to a 3d SDF file, out_f.
#     """
#     out_gen = (
#         StructGen.constrained_peptide_smiles(*pair)
#         for pair in process_constraints(fasta_f)
#     )
#     # pdb.set_trace()
#     StructGen.write_library(out_gen, out_f, write="structure", write_to_file=True)


# if __name__ == "__main__":
#     import argparse
#     parser = argparse.ArgumentParser(description="Generate peptides from a fasta file.")
#     parser.add_argument("fasta_f", help="Fasta file of peptides.")
#     parser.add_argument("--out_f", help="Output file.", default="out.sdf")
#     args = parser.parse_args()
    
#     main(args.fasta_f, args.out_f)
