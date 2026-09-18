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
from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import Crippen
import argparse


# Known synthesis difficulty patterns
forbidden_motifs = {
    "Over 2 prolines in a row are difficult to synthesise": r"[P]{3,}",
    "DG and DP are difficult to synthesise": r"D[GP]",
    "N or Q at N-terminus are difficult to synthesise": r"^[NQ]",
}

# List of charged residues
charged = ["H", "R", "K", "E", "D"]


@dataclass(frozen=True)
class SynthesisConfig:
    max_cysteines: int = 2
    max_glycine_run: int = 4
    max_length: int = 50
    charge_window: int = 5
    max_logp: float = 0.0


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


def check_glycine_runs(seq, max_run=4):
    # Check for runs of more than the allowed number of glycines
    pattern = rf"G{{{max_run + 1},}}"
    return (
        f"More than {max_run} consecutive glycines found"
        if re.search(pattern, seq)
        else None
    )


def check_length(seq, max_length=50):
    # Ensure sequence length does not exceed 50 residues
    return f"Peptide too long (length: {len(seq)})" if len(seq) > max_length else None


def check_charge(seq, window=5, charged_residues=None):
    # Check that there is at least one charged residue in every window-sized stretch
    charged_residues = charged if charged_residues is None else charged_residues
    if window < 1:
        raise ValueError("Charge window must be at least 1 residue")
    count = 0
    for resi in seq:
        count += 1
        if resi in charged_residues:
            count = 0
        if count >= window:
            return False
    return True


def collect_synthesis_issues(seq, config=None):
    # Aggregate all synthesis issues for a given sequence and SMILES
    config = SynthesisConfig() if config is None else config
    issues = check_forbidden_motifs(seq)
    # generate smiles from sequence
    try:
        # check if sequence can be converted to smiles
        smiles = Chem.MolToSmiles(Chem.MolFromSequence(seq))
        # Check hydrophobicity
        logp_val = log_partition_coefficient(smiles)
        if logp_val > config.max_logp:
            issues.append(
                f"Failed hydrophobicity: logP {logp_val:.2f} > {config.max_logp:.2f}"
            )
    except Exception:
        issues.append("Failed to generate SMILES, logP not checked.")

    # Check charge distribution
    if not check_charge(seq, window=config.charge_window):
        issues.append(
            f"Failed charge: need 1 charged residue every {config.charge_window} residues"
        )

    # Run additional structural checks
    for result in [
        check_cysteine_content(seq, max_c=config.max_cysteines),
        check_terminal_residues(seq),
        check_glycine_runs(seq, max_run=config.max_glycine_run),
        check_length(seq, max_length=config.max_length),
    ]:
        if result:
            issues.append(result)
    return issues


def evaluate_line(line, config=None):
    # Parse a single line (format: sequence-cyclization: smiles),
    # run all synthesis checks, and return issues or True if passed
    seq = line.strip()
    try:
        # pass if line contains '>'
        if ">" in line:
            return line, None

        # check if sequence has only natural amino acids
        if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
            raise ValueError("Line contains unnatural amino acids")

        issues = collect_synthesis_issues(seq, config=config)
        return (seq, True if not issues else issues)

    except Exception as e:
        # Return parsing error if line format is invalid
        return (seq or line, [f"Parsing error: {e}"])


def evaluate_file(input_file, output_file=None, config=None):
    # Evaluate all sequences; optionally write pass/fail results to a file
    results = []
    with open(input_file, "r") as f:
        for line in f:
            result = evaluate_line(line, config=config)
            results.append(result)

    # Write results to output file if provided
    if output_file is not None:
        header = ""
        with open(output_file, "w") as out:
            for line, result in results:
                if result is None:
                    # store lines starting with '>' as a variable to print later
                    header = line.strip()
                    continue
                if result is True:
                    out.write(f"{header}-synth_check_PASS\n{line.strip()}\n")
                else:
                    out.write(f"{header}-synth_check_FAIL:{result}\n{line.strip()}\n")
        return results

    # Print results to console
    else:
        header = ""
        for line, result in results:
            if result is None:
                header = line.strip()
                continue
            if result is True:
                print(f"{header} -> PASS\n{line.strip()}")
            else:
                print(f"{header} -> FAIL: {result}\n{line.strip()}")
        return results


def parse_args(argv=None):

    parser = argparse.ArgumentParser(
        description="Evaluate peptide synthesis feasibility from file."
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        help="Input FASTA file with natural peptide sequences.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=False,
        default=None,
        help="Optional output file to write results, otherwise output to terminal",
    )
    parser.add_argument(
        "--max_cysteines",
        type=int,
        default=2,
        help="Maximum allowed cysteines before the sequence is flagged.",
    )
    parser.add_argument(
        "--max_glycine_run",
        type=int,
        default=4,
        help="Maximum allowed consecutive glycines before the sequence is flagged.",
    )
    parser.add_argument(
        "--max_length",
        type=int,
        default=50,
        help="Maximum allowed peptide length.",
    )
    parser.add_argument(
        "--charge_window",
        type=int,
        default=5,
        help="Require at least one charged residue in every N-residue window.",
    )
    parser.add_argument(
        "--max_logp",
        type=float,
        default=0.0,
        help="Maximum allowed logP before the sequence is flagged as too hydrophobic.",
    )
    return parser.parse_args(argv)


def build_config(args):
    return SynthesisConfig(
        max_cysteines=args.max_cysteines,
        max_glycine_run=args.max_glycine_run,
        max_length=args.max_length,
        charge_window=args.charge_window,
        max_logp=args.max_logp,
    )


def main(argv=None):

    # CLI setup to specify input and output files
    parser = argparse.ArgumentParser(
        description="Evaluate peptide synthesis feasibility from file."
    )
    args = parse_args(argv)

    if args.max_cysteines < 0 or args.max_glycine_run < 0 or args.max_length < 0:
        parser.error("Threshold values must be non-negative")
    if args.charge_window < 1:
        parser.error("charge_window must be at least 1")

    # Run evaluation on the given file
    evaluate_file(args.input_file, args.output_file, config=build_config(args))


if __name__ == "__main__":
    main()
