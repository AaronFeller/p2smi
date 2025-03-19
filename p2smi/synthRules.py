import re

from rdkit import Chem
from rdkit.Chem import Crippen

forbidden_motifs = {
    "Over 2 prolines in a row are difficult to synthesise": r"[P]{3,}",
    "DG and DP are difficult to synthesise": r"D[GP]",
    "N or Q at N-terminus are difficult to synthesise": r"^[NQ]",
}
charged = ["H", "R", "K", "E", "D"]


class SmilesError(Exception):
    pass


def log_partition_coefficient(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SmilesError(f"Could not parse SMILES: {smiles}")
    return Crippen.MolLogP(mol)


def check_forbidden_motifs(seq):
    return [
        desc for desc, pattern in forbidden_motifs.items() if re.search(pattern, seq)
    ]


def check_cysteine_content(seq, max_c=2):
    c_count = seq.upper().count("C")
    return f"Too many cysteines (found {c_count})" if c_count > max_c else None


def check_terminal_residues(seq):
    return (
        "C-terminal residue is Pro or Cys, which may cause synthesis issues"
        if seq[-1] in ["P", "C"]
        else None
    )


def check_glycine_runs(seq):
    return (
        "More than 4 consecutive glycines found" if re.search(r"G{4,}", seq) else None
    )


def check_length(seq, max_length=50):
    return f"Peptide too long (length: {len(seq)})" if len(seq) > max_length else None


def check_charge(seq):
    count = 0
    for resi in seq:
        count += 1
        if resi in charged:
            count = 0
        if count >= 5:
            return False
    return True


def collect_synthesis_issues(smiles, seq):
    issues = check_forbidden_motifs(seq)
    try:
        logp_val = log_partition_coefficient(smiles)
        if logp_val > 0:
            issues.append(f"Failed hydrophobicity: logP {logp_val:.2f}")
    except SmilesError as e:
        issues.append(str(e))

    if not check_charge(seq):
        issues.append("Failed charge: need 1 charged residue every 5 residues")

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
    """Parse a line of format `sequence-cyclization: smiles`,
    run checks, return issues or True."""
    try:
        seq_part, smiles = line.strip().split(": ", 1)
        sequence, *_ = seq_part.split("-", 1)
        sequence = sequence.replace(",", "").strip().upper()
        issues = collect_synthesis_issues(smiles, sequence)
        return (line.strip(), True if not issues else issues)
    except Exception as e:
        return (line.strip(), [f"Parsing error: {e}"])


def evaluate_file(input_file, output_file=None):
    """Run synthesis checks on each line in a file
    and optionally write the results."""
    results = []
    with open(input_file, "r") as f:
        for line in f:
            result = evaluate_line(line)
            results.append(result)

    if output_file:
        with open(output_file, "w") as out:
            for line, result in results:
                if result is True:
                    out.write(f"{line} -> PASS\n")
                else:
                    out.write(f"{line} -> FAIL: {result}\n")
    print("Results:")
    for line, result in results:
        if result is True:
            print("PASS")
        else:
            print(f"FAIL: {result}")
    return results


if __name__ == "__main__":
    import argparse

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

    evaluate_file(args.input_file, args.output_file)
