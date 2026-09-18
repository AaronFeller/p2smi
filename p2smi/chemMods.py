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
from pathlib import Path
from rdkit import Chem
from rdkit import RDLogger

from p2smi.modifiers import (
    ModificationError,
    ModifierRegistry,
    SiteSelector,
    apply_modifier,
    apply_recipe,
    builtin_modifier_registry,
)

RDLogger.DisableLog("rdApp.*")  # quiet RDKit in batch

_AMIDE_N_PATTERN = Chem.MolFromSmarts(
    "[N;H1;X3;$([N]-[C](=O));$([N]-[C]-[C](=O))]"
)
_PEGYLATION_N_PATTERN = Chem.MolFromSmarts("[N;H1,H2;!$([N]-[C](=O))]")


def is_valid_smiles(smiles: str) -> bool:
    return Chem.MolFromSmiles(smiles) is not None


def _matching_atom_indices(mol, pattern):
    return sorted({match[0] for match in mol.GetSubstructMatches(pattern)})


def _pegylation_anchor_indices(mol):
    return _matching_atom_indices(mol, _PEGYLATION_N_PATTERN)


def add_n_methylation(sequence: str, methylation_residue_fraction: float):
    """
    Add methyl groups to peptide amide nitrogens selected via SMARTS matching.
    """
    if methylation_residue_fraction <= 0:
        return sequence, 0

    mol = Chem.MolFromSmiles(sequence)
    if mol is None:
        return sequence, 0

    match_indices = _matching_atom_indices(mol, _AMIDE_N_PATTERN)
    if not match_indices:
        return sequence, 0

    k = math.ceil(len(match_indices) * methylation_residue_fraction)
    chosen = random.sample(match_indices, min(k, len(match_indices)))

    product = mol
    definition = builtin_modifier_registry().modifiers["n_methyl"]
    for atom_idx in chosen:
        product = apply_modifier(
            product,
            definition,
            SiteSelector(atom_index=atom_idx),
        )

    return Chem.MolToSmiles(product, isomericSmiles=True), len(chosen)


def add_pegylation(sequence: str):
    """
    Acylate a random free amine with an mPEG chain containing 1-4 EO units.
    """
    mol = Chem.MolFromSmiles(sequence)
    if mol is None:
        return sequence, None

    anchors = _pegylation_anchor_indices(mol)
    if not anchors:
        return sequence, None

    anchor_idx = random.choice(anchors)
    peg_units = random.randint(1, 4)
    product = apply_modifier(
        mol,
        builtin_modifier_registry().modifiers["mpeg_acetyl"],
        SiteSelector(atom_index=anchor_idx),
        {"units": peg_units},
    )
    return Chem.MolToSmiles(product, isomericSmiles=True), "CCO" * peg_units


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


def load_modifier_registry(registry_file=None):
    registry = builtin_modifier_registry()
    if registry_file is not None:
        registry = ModifierRegistry.read_json(Path(registry_file), base=registry)
    return registry


def process_recipe_sequences(fp, recipe_id, registry, on_error="error"):
    if on_error not in {"error", "skip"}:
        raise ModificationError("on_error must be 'error' or 'skip'")
    try:
        recipe = registry.recipes[recipe_id]
    except KeyError as exc:
        raise ModificationError(f"Unknown recipe: {recipe_id}") from exc

    for header, sequence in parse_input_lines(fp):
        if sequence is None:
            if on_error == "error":
                raise ModificationError(f"Malformed input line: {header}")
            yield f"{header} [Skipped malformed line]"
            continue
        try:
            product, applied = apply_recipe(sequence, recipe, registry)
        except ModificationError as exc:
            if on_error == "error":
                raise
            yield f"{header or '[Unlabelled]'} [Skipped: {exc}]"
            continue
        annotation = f"[{' - '.join(applied)}]"
        yield f"{header}{annotation}: {product}" if header else f"{annotation}: {product}"


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


def process_recipe_file(input_file, output_file, recipe_id, registry, on_error="error"):
    with open(input_file, "r") as infile:
        lines = process_recipe_sequences(infile, recipe_id, registry, on_error=on_error)
        if output_file:
            with open(output_file, "w") as outfile:
                outfile.write("\n".join(lines))
        else:
            for line in lines:
                print(line)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Modify peptide SMILES randomly or with explicit recipes."
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
    parser.add_argument("--recipe", help="Deterministic recipe ID.")
    parser.add_argument(
        "--modifier-registry",
        type=Path,
        help="JSON modifier definitions and recipes to overlay on built-ins.",
    )
    parser.add_argument(
        "--on-error",
        choices=("error", "skip"),
        default="error",
        help="Recipe-mode behavior for malformed input or incompatible sites.",
    )
    args = parser.parse_args(argv)

    if args.recipe:
        try:
            process_recipe_file(
                args.input_file,
                args.output_file,
                args.recipe,
                load_modifier_registry(args.modifier_registry),
                on_error=args.on_error,
            )
        except ModificationError as exc:
            parser.error(str(exc))
        return

    process_file(
        args.input_file,
        args.output_file,
        args.peg_rate,
        args.nmeth_rate,
        args.nmeth_residues,
    )


if __name__ == "__main__":
    main()
