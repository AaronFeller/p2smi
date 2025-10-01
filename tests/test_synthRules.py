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
    issues = collect_synthesis_issues("FWYLIV")
    assert any("logP" in issue for issue in issues)


def test_failure_to_generate_smiles():
    issues = collect_synthesis_issues("XXXX")
    assert any(
        "Failed to generate SMILES, logP not checked." in issue for issue in issues
    )


def test_evaluate_header_pass_case():
    line = ">seq1"
    result = evaluate_line(line)
    assert result[0] == ">seq1"
    assert result[1] is None


def test_evaluate_line_pass_case():
    line = "AARIN"
    result = evaluate_line(line)
    assert result[0] == "AARIN"
    assert result[1] is True


def test_evaluate_line_fail_case_forbidden_motif():
    line = "PPP"
    result = evaluate_line(line)
    assert any("prolines" in issue.lower() for issue in result[1])


def test_evaluate_file(tmp_path):
    test_content = ">Test1\n" "ASK\n" ">Test2\n" "PPP\n"
    test_file = tmp_path / "test_input.txt"
    output_file = tmp_path / "test_output.txt"
    test_file.write_text(test_content)

    results = evaluate_file(test_file, output_file)
    assert len(results) == 4
    assert results[0][1] is None
    assert results[1][1] is True
    assert results[2][1] is None
    assert type(results[3][1]) is list

    written_content = output_file.read_text()
    assert "PASS" in written_content
    assert "FAIL" in written_content
