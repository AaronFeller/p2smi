from p2smi.synthRules import (
    check_charge,
    check_cysteine_content,
    check_forbidden_motifs,
    collect_synthesis_issues,
    evaluate_file,
    evaluate_line,
)


def test_check_forbidden_motifs_detects_patterns():
    assert (
        "Over 2 prolines in a row are difficult to synthesise"
        in check_forbidden_motifs("APPPG")
    )
    assert "DG and DP are difficult to synthesise" in check_forbidden_motifs("DPQ")
    assert "N or Q at N-terminus are difficult to synthesise" in check_forbidden_motifs(
        "NAGK"
    )


def test_check_cysteine_content_limit():
    assert check_cysteine_content("CCCC") == "Too many cysteines (found 4)"
    assert check_cysteine_content("CC") is None


def test_test_charge_behavior():
    assert check_charge("KRRDK") is True
    assert check_charge("AAAAAAH") is False  # no charged residue in >5 stretch


def test_collect_synthesis_issues_catches_hydrophobicity():
    # Very hydrophobic simple molecule
    smiles = "CCCCCCCCCCCC"
    issues = collect_synthesis_issues(smiles, "ACDEFG")
    assert any("logP" in issue for issue in issues)


def test_evaluate_line_pass_case():
    line = "A,S,K-LINEAR: N[C@@H](C)C(=O)N[C@@H]" "(CCCCN)C(=O)N[C@@H](C(C)C)C(=O)O"
    result = evaluate_line(line)
    assert result[1] is True


def test_evaluate_line_fail_case_forbidden_motif():
    line = "P,P,P-LINEAR: N[C@@H](C)C(=O)N[C@@H]" "(C(C)C)C(=O)N[C@@H](C(C)C)C(=O)O"
    result = evaluate_line(line)
    assert any("prolines" in issue.lower() for issue in result[1])


def test_evaluate_file(tmp_path):
    test_content = (
        "A,S,K-LINEAR: N[C@@H](C)C(=O)N[C@@H]"
        "(CCCCN)C(=O)N[C@@H](C(C)C)C(=O)O\n"
        "P,P,P-LINEAR: N[C@@H](C)C(=O)N[C@@H]"
        "(C(C)C)C(=O)N[C@@H](C(C)C)C(=O)O"
    )
    test_file = tmp_path / "test_input.txt"
    output_file = tmp_path / "test_output.txt"
    test_file.write_text(test_content)

    results = evaluate_file(test_file, output_file)
    assert len(results) == 2
    assert results[0][1] is True
    assert any("prolines" in issue.lower() for issue in results[1][1])

    written_content = output_file.read_text()
    assert "PASS" in written_content
    assert "FAIL" in written_content
