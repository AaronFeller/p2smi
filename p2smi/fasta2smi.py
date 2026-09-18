"""
Module to input FASTA peptide files and generate 3D structures.

Original by Fergal; modified by Aaron Feller (2025).
Reads peptide sequences from a FASTA file, resolves structural constraints,
converts them to SMILES strings, and writes outputs.

Uses p2smi.utilities.smilesgen.
"""

import argparse
import p2smi.utilities.smilesgen as smilesgen
from p2smi.modifiers import (
    ModificationError,
    ModifierRegistry,
    build_modified_peptide,
    builtin_modifier_registry,
)
from p2smi.registry import RegistryError, ResidueRegistry, parse_sequence_tokens

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


def constraint_resolver(sequence, constraint, registry=None):
    # Resolve constraints by checking known patterns or fallback attempts.
    # Return (sequence, constraint) or fallback to linear if none apply.
    functions = {
        "SS": smilesgen.can_ssbond,
        "HT": smilesgen.can_htbond,
        "SCNT": smilesgen.can_scntbond,
        "SCCT": smilesgen.can_scctbond,
        "SCSC": smilesgen.can_scscbond,
    }
    constraint_functions = {
        key: (
            (lambda seq, func=func: func(seq))
            if registry is None
            else (lambda seq, func=func: func(seq, registry=registry))
        )
        for key, func in functions.items()
    }

    if constraint.upper() in constraint_functions:
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


def validate_constraint_pattern(peptideseq, pattern, registry=None):
    """
    Validate that a user-supplied constraint pattern matches both
    the peptide chemistry (using all_aminos) and the positional mask.

    Accepts 'N'/'E' sidechains that link to the peptide C-terminus and 'Z'
    sidechains that link to the peptide N-terminus.
    Returns (is_valid, message) instead of raising exceptions.
    """

    if registry is not None:
        try:
            seq = parse_sequence_tokens(peptideseq, registry)
            residues = [registry.resolve(token) for token in seq]
            capabilities = []
            for residue in residues:
                mol = smilesgen._mol_from_smiles(
                    residue.smiles, context=f"residue {residue.id}"
                )
                sites = smilesgen._reactive_site_indices(mol)
                capabilities.append(
                    {kind for kind, indices in sites.items() if len(indices) == 1}
                )
        except (RegistryError, smilesgen.SmilesError) as exc:
            return False, str(exc)
    else:
        seq = "".join(peptideseq) if not isinstance(peptideseq, str) else peptideseq
        try:
            legacy_residues = [all_aminos[LETTER2NAME[r]] for r in seq]
        except KeyError as exc:
            return False, f"Undefined residue {exc.args[0]} in sequence."
        residues = legacy_residues
        capabilities = [
            {key for key in ("disulphide", "nterm", "ester", "cterm") if has_capability(aa, key)}
            for aa in residues
        ]
    tag = pattern[:2].upper()
    mask = pattern[2:]

    if not pattern:
        return True, "No constraint pattern provided."

    if len(mask) != len(seq) and len(mask) > 0:
        return (
            False,
            f"Mask length ({len(mask)}) does not match sequence length ({len(seq)}).",
        )

    # -----------------------------
    # Validate per-position codes
    # -----------------------------
    for i, code in enumerate(mask):
        if code == "X":
            continue
        residue_capabilities = capabilities[i]
        if code == "C" and "disulphide" not in residue_capabilities:
            return False, f"Position {i}: '{seq[i]}' cannot form disulfide ('C')."
        if code == "N" and "nterm" not in residue_capabilities:
            return False, f"Position {i}: '{seq[i]}' lacks nterm capability ('N')."
        if code == "E" and "ester" not in residue_capabilities:
            return False, f"Position {i}: '{seq[i]}' lacks ester capability ('E')."
        if code == "Z" and "cterm" not in residue_capabilities:
            return False, f"Position {i}: '{seq[i]}' lacks cterm capability ('Z')."
        if code not in {"X", "C", "N", "E", "Z"}:
            return False, f"Position {i}: unknown constraint code '{code}'."

    # -----------------------------
    # Constraint-type validation
    # -----------------------------
    if tag == "SS":
        if mask.count("C") != 2 or set(mask) - {"X", "C"}:
            return False, "SS constraint mask must select exactly two disulphide sites."
        return True, "Valid SS constraint."

    elif tag == "HT":
        if mask:
            return False, "HT constraint does not accept a positional mask."
        if len(residues) < 2:
            return False, "HT requires at least two residues."
        return True, "Valid HT constraint."

    elif tag == "SC":
        codes = set(mask) - {"X"}
        if mask.count("Z") > 1 or mask.count("N") + mask.count("E") > 1:
            return False, "SC constraint must select at most one acid and one donor site."
        if ({"N", "E"} & codes) and "Z" in codes:
            subtype = "SCSC"  # sidechain–sidechain
        elif {"N", "E"} & codes:
            subtype = "SCCT"  # sidechain–C-terminus
        elif "Z" in codes:
            subtype = "SCNT"  # sidechain–N-terminus
        else:
            return False, f"Unrecognized SC pattern: {pattern}"

        # subtype-specific residue chemistry beyond the per-position checks above
        if subtype == "SCSC":
            n_like = any(code in {"N", "E"} for code in mask)
            c_like = "Z" in codes
            if not (n_like and c_like):
                return False, "SCSC requires one N/E site and one Z site in the mask."

        return True, f"Valid {subtype} constraint."

    else:
        return False, f"Unknown constraint tag '{tag}'."


def normalize_constraint(result):
    # if it's a tuple like (seq, pattern)
    if isinstance(result, tuple):
        return result[1]
    return result


def process_constraints(fasta_file, registry=None):
    return (
        (
            seq,
            (
                constr
                if "X" in constr
                else normalize_constraint(constraint_resolver(seq, constr, registry=registry))
            ),
        )
        for seq, constr in parse_fasta(fasta_file)
    )


def generate_smiles_strings(
    input_fasta,
    out_file,
    verbose=False,
    registry_file=None,
    modifier_registry_file=None,
    recipe_id=None,
    on_error="error",
):
    registry = None
    if registry_file is not None:
        registry = ResidueRegistry.from_legacy_mapping(all_aminos).merge(
            ResidueRegistry.read_json(registry_file)
        )
    resolved_sequences = list(process_constraints(input_fasta, registry=registry))
    valid_sequences = []

    for seq, constr in resolved_sequences:
        is_valid, message = validate_constraint_pattern(seq, constr, registry=registry)
        if verbose:
            print(f"[DEBUG] {seq=} {constr=} -> {is_valid=} {message=}")
        if not is_valid:
            print(f"Warning: {message}")
            continue
        valid_sequences.append((seq, constr))

    if recipe_id:
        residue_registry = registry or ResidueRegistry.from_legacy_mapping(all_aminos)
        modifier_registry = builtin_modifier_registry()
        if modifier_registry_file:
            modifier_registry = ModifierRegistry.read_json(
                modifier_registry_file, base=modifier_registry
            )
        try:
            recipe = modifier_registry.recipes[recipe_id]
        except KeyError as exc:
            raise ModificationError(f"Unknown recipe: {recipe_id}") from exc

        def modified_entries():
            for seq, constr in valid_sequences:
                try:
                    product, applied = build_modified_peptide(
                        seq,
                        constr,
                        residue_registry,
                        modifier_registry,
                        recipe,
                    )
                except Exception as exc:
                    if on_error == "error":
                        raise
                    print(f"Warning: skipped {seq}: {exc}")
                    continue
                bond_label = f"{constr or 'linear'}[{','.join(applied)}]"
                yield seq, bond_label, product

        entries = modified_entries()
    else:
        entries = (
            (
                smilesgen.constrained_peptide_smiles(seq, constr)
                if registry is None
                else smilesgen.constrained_peptide_smiles(seq, constr, registry=registry)
            )
            for seq, constr in valid_sequences
        )
    smilesgen.write_library(
        entries,
        out_file,
        write="text",
        write_to_file=True,
    )


def main(argv=None):
    # CLI entry point: takes FASTA file input, output file path, and generates structures
    parser = argparse.ArgumentParser(
        description=(
            "Convert peptide FASTA files into SMILES strings with optional structural "
            "constraints.\n\n"
            "Each FASTA entry should use the notation:\n"
            "  >peptide_name|CONSTRAINT\n"
            "  ACDEFGHIKLMNPQRSTVWY\n\n"
            "  If CONSTRAINT is not used, peptide will be treated as linear.\n\n"
            "  To let the program infer residues for CONSTRAINT, use:\n"
            "    '>PEPTIDE_NAME|{SS,HT,SCSC,SCNT,SCCT}'.\n\n"
            "  To define cyclization residues manually, encode pattern using:\n"
            "    X    - Any residue\n"
            "    C    - Cysteine (for disulfide bonds)\n"
            "    N    - Nucleophilic sidechain bonded to C-terminal (e.g., K, S, T, Y, C)\n"
            "    E    - Ester sidechain bonded to C-terminal (optional advanced form)\n"
            "    Z    - Carboxyl sidechain bonded to N-terminal (e.g., D, E)\n\n"
            "  Example manual CONSTRAINT:\n"
            "    SS   - SSXXXXCXXXCX (Disulfide bond)\n"
            "    SCSC - SCXXNXXXXZ (Sidechain–sidechain linkage)\n"
            "    SCNT - SCXXXXXZXXX (Sidechain–N-terminus linkage)\n"
            "    SCCT - SCXXNXXXXXX (Sidechain–C-terminus linkage)\n\n"
            "Residue capabilities are validated using the amino acid database"
            " in p2smi.utilities.smilesgen.\n"
            "If an invalid or incompatible pattern is detected, a warning is"
            " printed and the peptide is skipped.\n\n"
            "For documentation and examples, visit: https://github.com/aaronfeller/p2smi"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input_fasta", required=True, help="FASTA file of peptides."
    )
    parser.add_argument("-o", "--out_file", required=True, help="Output file.")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose output."
    )
    parser.add_argument(
        "--registry-file",
        help="Optional v2 JSON residue registry for comma-delimited sequence IDs.",
    )
    parser.add_argument(
        "--recipe",
        help="Explicit modification recipe ID to apply during sequence assembly.",
    )
    parser.add_argument(
        "--modifier-registry",
        help="JSON modifier definitions and recipes to overlay on built-ins.",
    )
    parser.add_argument(
        "--on-error",
        choices=("error", "skip"),
        default="error",
        help="Recipe-mode behavior for incompatible modification sites.",
    )
    args = parser.parse_args(argv)

    generate_smiles_strings(
        args.input_fasta,
        args.out_file,
        args.verbose,
        registry_file=args.registry_file,
        modifier_registry_file=args.modifier_registry,
        recipe_id=args.recipe,
        on_error=args.on_error,
    )


if __name__ == "__main__":
    main()
