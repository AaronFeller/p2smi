from p2smi.chemMods import (
    add_n_methylation,
    add_pegylation,
    is_valid_smiles,
    modify_sequence,
    parse_input_lines,
    process_sequences,
)


def test_is_valid_smiles():
    assert is_valid_smiles("C1CCCCC1")  # cyclohexane
    assert not is_valid_smiles("INVALID_SMILES")


def test_add_n_methylation_counts():
    test_seq = "C(=O)N[C@H]C(=O)N[C@H]C"  # 2 peptidic bonds
    modified_seq, count = add_n_methylation(test_seq, 1.0)
    assert count == 2
    assert modified_seq.count("(C)") == 2


def test_add_pegylation_adds_peg():
    test_seq = "N[C@@H](CCCCN)C(=O)O"  # L-Lysine
    peg_seq, peg = add_pegylation(test_seq)
    assert peg_seq != test_seq
    assert peg and peg.startswith("O") and "CCO" in peg


def test_parse_input_lines_handles_good_and_bad_lines():
    lines = [
        "peptide1: C1CCCCC1",
        "malformed line without colon",
        "peptide2: C1CCNCC1",
    ]
    parsed = list(parse_input_lines(lines))
    assert parsed[0] == ("peptide1", "C1CCCCC1")
    assert parsed[1] == ("[Malformed line]", None)
    assert parsed[2] == ("peptide2", "C1CCNCC1")


def test_modify_sequence_modifications_included():
    test_seq = "C(=O)N[C@H]C(=O)N[C@H]C"
    mod_seq, mods = modify_sequence(
        test_seq, methylate=True, pegylate=True, nmeth_residues=0.5
    )
    assert any("N-methylation" in m for m in mods)
    assert any("PEGylation" in m for m in mods)


def test_process_sequences_runs_and_marks_modifications():
    input_lines = [
        "pep1: C(=O)N[C@H]C(=O)N[C@H]C",
        "pep2: N[C@@H](CCCCN)C(=O)O",
        "bad_line",
    ]
    results = list(
        process_sequences(input_lines, nmeth_rate=1.0, peg_rate=1.0, nmeth_residues=0.5)
    )
    assert any("pep1" in r and "N-methylation" in r for r in results)
    assert any("pep2" in r and "PEGylation" in r for r in results)
    assert any("Skipped malformed line" in r for r in results)
