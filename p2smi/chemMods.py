#!/usr/bin/env python
"""
Fast peptide SMILES modifier (PEGylation + N-methylation), streamed.

Speedups:
- Stream input lines; no full-file materialization.
- Bernoulli per-line decisions (no preselected index sets).
- Precompiled regex patterns reused across lines.
- O(L + #inserts) builders for string edits.
- Single RDKit validation per final sequence.
"""

import argparse
import math
import random
import re
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")  # quiet RDKit in batch

# --- Precompiled patterns ---
# Amide N to N-methylate: C(=O)N[C@   (insert "(C)" after the 'N')
_AMIDE_N_PATTERN = re.compile(r"C\(=O\)N\[C@")
# PEG insertion anchor: ... "CN)"  (insert PEG after 'CN')
_PEG_ANCHOR_PATTERN = re.compile(r"CN\)")


def is_valid_smiles(smiles: str) -> bool:
    return Chem.MolFromSmiles(smiles) is not None


def _insert_many(base: str, inserts):
    """
    Insert multiple (idx, text) into base in a single pass.
    `inserts` must be sorted by idx ascending.
    """
    if not inserts:
        return base
    out = []
    prev = 0
    for idx, txt in inserts:
        out.append(base[prev:idx])
        out.append(txt)
        prev = idx
    out.append(base[prev:])
    return "".join(out)


def add_n_methylation(sequence: str, methylation_residue_fraction: float):
    """
    Insert '(C)' after the amide N for a random subset of matches.
    We compute insertion indices once and build result in one pass.
    """
    if methylation_residue_fraction <= 0:
        return sequence, 0

    # Find starts of the pattern; insert after the 'N' (offset +6)
    # 'C(=O)N' is 6 chars before the '['
    match_starts = [m.start() for m in _AMIDE_N_PATTERN.finditer(sequence)]
    if not match_starts:
        return sequence, 0

    k = math.ceil(len(match_starts) * methylation_residue_fraction)
    chosen = random.sample(match_starts, min(k, len(match_starts)))
    # Build list of (absolute) insertion points
    inserts = sorted(((pos + 6, "(C)") for pos in chosen), key=lambda x: x[0])

    return _insert_many(sequence, inserts), len(chosen)


def add_pegylation(sequence: str):
    """
    Insert a random-length PEG chain after a random 'CN)' anchor.
    PEG = O(CCO){1..4}C  (length picked uniformly)
    """
    anchors = [m.start() for m in _PEG_ANCHOR_PATTERN.finditer(sequence)]
    if not anchors:
        return sequence, None

    pos = random.choice(anchors)
    peg = "O" + "CCO" * random.randint(1, 4) + "C"
    # Insert right after 'CN' (i.e., after pos+2)
    insert_idx = pos + 2
    return _insert_many(sequence, [(insert_idx, peg)]), peg


def parse_input_lines(fp):
    """
    Stream parser: yields (header, smiles) or ("[Malformed line]", None).
    Accepts either 'Header: SMILES' or just 'SMILES' per line.
    """
    for raw in fp:
        line = raw.strip()
        if not line:
            continue
        if ": " in line:
            header, smi = line.split(": ", 1)
            header, smi = header.strip(), smi.strip()
            # don’t validate here; do it once after modification
            yield (header, smi)
        else:
            # bare SMILES or malformed
            yield (
                ("", line)
                if Chem.MolFromSmiles(line) is not None
                else ("[Malformed line]", None)
            )


def modify_sequence(
    sequence: str, do_methylate: bool, do_pegylate: bool, nmeth_residues: float
):
    mods = []
    seq = sequence

    if do_methylate:
        seq, methyl_count = add_n_methylation(seq, nmeth_residues)
        mods.append(f"N-methylation({methyl_count})")

    if do_pegylate:
        seq, peg = add_pegylation(seq)
        if peg:
            mods.append(f"PEGylation({peg.count('CCO')})")
        else:
            mods.append("PEGylation('N/A')")

    return seq, mods


def process_sequences(fp, nmeth_rate: float, peg_rate: float, nmeth_residues: float):
    """
    Stream through file-like fp; decide per line via Bernoulli(p),
    apply edits, validate once, and yield output strings.
    """
    for header, seq in parse_input_lines(fp):
        if seq is None:
            yield f"{header} [Skipped malformed line]"
            continue

        # Bernoulli per line (fast; no preselect sets)
        do_methylate = (random.random() < nmeth_rate) if nmeth_rate > 0 else False
        do_pegylate = (random.random() < peg_rate) if peg_rate > 0 else False

        mod_seq, mods = modify_sequence(seq, do_methylate, do_pegylate, nmeth_residues)

        if not is_valid_smiles(mod_seq):
            yield f"{header} [Invalid SMILES skipped]"
            continue

        mod_str = f"[{' - '.join(mods)}]" if mods else ""
        prefix = f"{header}{mod_str}".strip()
        # Preserve your original "Header: SMILES" shape; allow empty header
        if header:
            yield f"{prefix}: {mod_seq}"
        else:
            yield f"{mod_seq}" if not mods else f"{mod_str}: {mod_seq}"


def process_file(
    input_file: str,
    output_file: str,
    peg_rate: float,
    nmeth_rate: float,
    nmeth_residues: float,
):
    if output_file:
        with open(input_file, "r") as infile, open(output_file, "w") as outfile:
            first = True
            for line in process_sequences(infile, nmeth_rate, peg_rate, nmeth_residues):
                if not first:
                    outfile.write("\n")
                outfile.write(line)
                first = False
            if first:
                outfile.write("")  # no lines
    else:
        with open(input_file, "r") as infile:
            for line in process_sequences(infile, nmeth_rate, peg_rate, nmeth_residues):
                print(line)


def main():
    parser = argparse.ArgumentParser(
        description="Modify peptide SMILES with PEGylation and N-methylation."
    )
    parser.add_argument("-i", "--input_file", required=True, help="Input file path.")
    parser.add_argument(
        "-o",
        "--output_file",
        help="Optional output file path; prints to stdout if omitted.",
    )
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
        help="Fraction of amide sites per sequence to N-methylate (0-1).",
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
