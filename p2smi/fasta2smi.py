"""
Module to input FASTA peptide files and generate 3D structures.

Original by Fergal; modified by Aaron Feller (2025).
Reads peptide sequences from a FASTA file, resolves structural constraints,
converts them to SMILES strings, and writes outputs.

Uses p2smi.utilities.smilesgen.
"""

import argparse
import p2smi.utilities.smilesgen as smilesgen
all_aminos = smilesgen.all_aminos
LETTER2NAME = smilesgen.LETTER2NAME

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

def has_capability(aa, key):
    """
    Return True if the amino acid supports the given capability key.
    Works whether the value is a bool, 'False' string, or a SMILES fragment.
    """
    val = aa.get(key, False)
    if val is None:
        return False
    if isinstance(val, str):
        # strip whitespace and quotes, and check for the literal word False
        val_stripped = val.strip().lower()
        if val_stripped in {"", "false", "none"}:
            return False
        return True
    return bool(val)

def validate_constraint_pattern(peptideseq, pattern):
    """
    Validate that a user-supplied constraint pattern matches both
    the peptide chemistry (using all_aminos) and the positional mask.

    Uses only 'N' for nucleophilic sidechains (no 'E').
    Returns (is_valid, message) instead of raising exceptions.
    """

    seq = "".join(peptideseq) if not isinstance(peptideseq, str) else peptideseq
    tag = pattern[:2].upper()
    mask = pattern[2:]

    if not pattern:
        return True, "No constraint pattern provided."

    if len(mask) != len(seq) and len(mask) > 0:
        return False, f"Mask length ({len(mask)}) does not match sequence length ({len(seq)})."

    try:
        residues = [all_aminos[LETTER2NAME[r]] for r in seq]
    except KeyError as e:
        return False, f"Undefined residue {e.args[0]} in sequence."

    # -----------------------------
    # Validate per-position codes
    # -----------------------------
    for i, code in enumerate(mask):
        if code == "X":
            continue
        aa = residues[i]
        if code == "C" and not has_capability(aa, "disulphide"):
            return False, f"Position {i}: '{seq[i]}' cannot form disulfide ('C')."
        if code == "N" and not has_capability(aa, "nterm"):
            return False, f"Position {i}: '{seq[i]}' lacks nterm capability ('N')."
        if code == "Z" and not has_capability(aa, "cterm"):
            return False, f"Position {i}: '{seq[i]}' lacks cterm capability ('Z')."

    # -----------------------------
    # Constraint-type validation
    # -----------------------------
    if tag == "SS":
        if sum(has_capability(aa, "disulphide") for aa in residues) < 2:
            return False, "SS constraint requires ≥2 disulphide-capable residues."
        return True, "Valid SS constraint."

    elif tag == "HT":
        if len(residues) < 2:
            return False, "HT requires at least two residues."
        return True, "Valid HT constraint."

    elif tag == "SC":
        codes = set(mask)
        if "N" in codes and "Z" in codes:
            subtype = "SCSC"  # sidechain–sidechain
        elif "N" in codes:
            subtype = "SCNT"  # sidechain–N-terminus (now covers previous 'E')
        elif "Z" in codes:
            subtype = "SCCT"  # sidechain–C-terminus
        else:
            return False, f"Unrecognized SC pattern: {pattern}"

        # subtype-specific residue chemistry
        if subtype == "SCSC":
            n_like = any(has_capability(aa, "nterm") for aa in residues)
            c_like = any(has_capability(aa, "cterm") for aa in residues)
            if not (n_like and c_like):
                return False, "SCSC requires one nterm=True and one cterm=True residue."
        elif subtype == "SCNT":
            if not any(has_capability(aa, "nterm") for aa in residues):
                return False, "SCNT requires at least one nterm=True residue."
        elif subtype == "SCCT":
            if not any(has_capability(aa, "cterm") for aa in residues):
                return False, "SCCT requires at least one cterm=True residue."

        return True, f"Valid {subtype} constraint."

    else:
        return False, f"Unknown constraint tag '{tag}'."

def normalize_constraint(result):
    # if it's a tuple like (seq, pattern)
    if isinstance(result, tuple):
        return result[1]
    return result

def process_constraints(fasta_file):
    return (
        (
            seq, constr if "X" in constr else normalize_constraint(
                constraint_resolver(seq, constr)
            )
        )
        for seq, constr in parse_fasta(fasta_file)
    )

def generate_smiles_strings(input_fasta, out_file, verbose=False):
    resolved_sequences = list(process_constraints(input_fasta))  # <-- materialize

    for seq, constr in resolved_sequences:
        is_valid, message = validate_constraint_pattern(seq, constr)
        if verbose:
            print(f"[DEBUG] {seq=} {constr=} -> {is_valid=} {message=}")
        if not is_valid:
            print(f"Warning: {message}")
            continue

    smilesgen.write_library(
        (smilesgen.constrained_peptide_smiles(seq, constr)
        for seq, constr in resolved_sequences),
        out_file,
        write="text",
        write_to_file=True,
    )


def main():
    # CLI entry point: takes FASTA file input, output file path, and generates structures
    parser = argparse.ArgumentParser(
        description=(
            "Convert peptide FASTA files into SMILES strings with optional structural constraints.\n\n"
            "Each FASTA entry should use the notation:\n"
            "  >peptide_name|CONSTRAINT\n"
            "  ACDEFGHIKLMNPQRSTVWY\n\n"
            
            "  If CONSTRAINT is not used, peptide will be treated as linear.\n\n"
            "  To let the program infer residues for CONSTRAINT, use:\n" 
            "    '>PEPTIDE_NAME|{SS,HT,SCSC,SCNT,SCCT}'.\n\n"
            "  To define cyclization residues manually, encode pattern using:\n"
            "    X    - Any residue\n"
            "    C    - Cysteine (for disulfide bonds)\n"
            "    N    - Residue bonded to N-terminal (e.g., K, S, T, Y, C)\n"
            "    Z    - Residue bonded to C-terminal (e.g., D, E, K, R)\n\n" 
            
            "  Example manual CONSTRAINT:\n"
            "    SS   - SSXXXXCXXXCX (Disulfide bond)\n"
            "    SCSC - SCXXNXXXXZ (Sidechain–sidechain linkage)\n"
            "    SCNT - SCXXNXXXXXX (Sidechain–N-terminus linkage)\n"
            "    SCCT - SCXXXXXZXXX (Sidechain–C-terminus linkage)\n\n"
            "Residue capabilities are validated using the amino acid database in p2smi.utilities.smilesgen.\n"
            "If an invalid or incompatible pattern is detected, a warning is printed and the peptide is skipped.\n\n"
            "For documentation and examples, visit: https://github.com/aaronfeller/p2smi"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
)
    parser.add_argument("-i", "--input_fasta", required=True, help="FASTA file of peptides.")
    parser.add_argument("-o", "--out_file", required=True, help="Output file.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output.")
    args = parser.parse_args()

    generate_smiles_strings(args.input_fasta, args.out_file, args.verbose)


if __name__ == "__main__":
    main()
