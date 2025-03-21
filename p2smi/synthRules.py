"""
Peptide Synthesis Feasibility Evaluator

This script reads peptide sequences (with associated SMILES) and evaluates their
synthetic feasibility based on a variety of rules:
- Detects forbidden motifs (e.g. proline runs, DG or DP sequences, N/Q at N-terminus)
- Checks cysteine content (too many cysteines can complicate synthesis)
- Warns if terminal residues are Pro or Cys
- Checks for long glycine runs (over 4 in a row)
- Ensures peptide length does not exceed recommended limits
- Validates hydrophobicity (logP should not be overly high)
- Checks charge distribution (at least one charged residue every 5 residues)

Input format (per line):
`sequence-cyclization: smiles`
(for an example input format, run genPeps.py)

Outputs pass/fail results along with reasoned diagnostics.
Uses RDKit for chemical property calculations.
"""

import re
from rdkit import Chem
from rdkit.Chem import Crippen

# Known synthesis difficulty patterns
forbidden_motifs = {
    "Over 2 prolines in a row are difficult to synthesise": r"[P]{3,}",
    "DG and DP are difficult to synthesise": r"D[GP]",
    "N or Q at N-terminus are difficult to synthesise": r"^[NQ]",
}

# List of charged residues
charged = ["H", "R", "K", "E", "D"]


class SmilesError(Exception):
    # Custom error for invalid SMILES strings
    pass


def log_partition_coefficient(smiles):
    # Calculate logP from a SMILES string; raise an error if invalid
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SmilesError(f"Could not parse SMILES: {smiles}")
    return Crippen.MolLogP(mol)


def check_forbidden_motifs(seq):
    # Check the sequence for forbidden motifs and return descriptive messages
    return [
        desc for desc, pattern in forbidden_motifs.items() if re.search(pattern, seq)
    ]


def check_cysteine_content(seq, max_c=2):
    # Check if cysteine count exceeds the allowed maximum
    c_count = seq.upper().count("C")
    return f"Too many cysteines (found {c_count})" if c_count > max_c else None


def check_terminal_residues(seq):
    # Warn if the C-terminal residue is Proline or Cysteine
    return (
        "C-terminal residue is Pro or Cys, which may cause synthesis issues"
        if seq[-1] in ["P", "C"]
        else None
    )


def check_glycine_runs(seq):
    # Check for runs of more than 4 glycines
    return (
        "More than 4 consecutive glycines found" if re.search(r"G{4,}", seq) else None
    )


def check_length(seq, max_length=50):
    # Ensure sequence length does not exceed 50 residues
    return f"Peptide too long (length: {len(seq)})" if len(seq) > max_length else None


def check_charge(seq):
    # Check that there is at least one charged residue every 5 residues
    count = 0
    for resi in seq:
        count += 1
        if resi in charged:
            count = 0
        if count >= 5:
            return False
    return True


def collect_synthesis_issues(smiles, seq):
    # Aggregate all synthesis issues for a given sequence and SMILES
    issues = check_forbidden_motifs(seq)
    try:
        # Check hydrophobicity
        logp_val = log_partition_coefficient(smiles)
        if logp_val > 0:
            issues.append(f"Failed hydrophobicity: logP {logp_val:.2f}")
    except SmilesError as e:
        issues.append(str(e))

    # Check charge distribution
    if not check_charge(seq):
        issues.append("Failed charge: need 1 charged residue every 5 residues")

    # Run additional structural checks
    for check_fn in [
        check_cysteine_content,
        check_terminal_residues,
        check_glycine_runs,
        check_length,
    ]:
        result = check_fn(seq)
        if result:
            issues.append(result)
    return issues


def evaluate_line(line):
    # Parse a single line (format: sequence-cyclization: smiles),
    # run all synthesis checks, and return issues or True if passed
    try:
        seq_part, smiles = line.strip().split(": ", 1)
        sequence, *_ = seq_part.split("-", 1)
        sequence = sequence.replace(",", "").strip().upper()
        issues = collect_synthesis_issues(smiles, sequence)
        return (line.strip(), True if not issues else issues)
    except Exception as e:
        # Return parsing error if line format is invalid
        return (line.strip(), [f"Parsing error: {e}"])


def evaluate_file(input_file, output_file=None):
    # Evaluate all sequences; optionally write pass/fail results to a file
    results = []
    with open(input_file, "r") as f:
        for line in f:
            result = evaluate_line(line)
            results.append(result)

    # Write results to output file if provided
    if output_file:
        with open(output_file, "w") as out:
            for line, result in results:
                if result is True:
                    out.write(f"{line} -> PASS\n")
                else:
                    out.write(f"{line} -> FAIL: {result}\n")

    # Print results to console
    print("Results:")
    for line, result in results:
        if result is True:
            print("PASS")
        else:
            print(f"FAIL: {result}")
    return results


def main():
    import argparse

    # CLI setup to specify input and output files
    parser = argparse.ArgumentParser(
        description="Evaluate peptide synthesis feasibility from file."
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        help="Input file with peptide-smiles lines",
    )
    parser.add_argument(
        "-o", "--output_file", help="Optional output file to write results"
    )
    args = parser.parse_args()

    # Run evaluation on the given file
    evaluate_file(args.input_file, args.output_file)


if __name__ == "__main__":
    main()
