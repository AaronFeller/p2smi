import pytest

# Mock smilesgen
import p2smi.utilities.smilesgen as smilesgen
from p2smi.fasta2smi import (
    InvalidConstraintError,
    constraint_resolver,
    parse_fasta,
    process_constraints,
)

smilesgen.can_ssbond = lambda seq: (seq, "SS")
smilesgen.can_htbond = lambda seq: (seq, "HT")
smilesgen.can_scntbond = lambda seq: (seq, "SCNT")
smilesgen.can_scctbond = lambda seq: (seq, "SCCT")
smilesgen.can_scscbond = lambda seq: (seq, "SCSC")
smilesgen.what_constraints = lambda seq: ["SS", "HT"]


def test_parse_fasta(tmp_path):
    fasta_content = ">seq1|SS\nACDE\n>seq2\nFGHI"
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    results = list(parse_fasta(fasta_file))
    assert results == [("ACDE", "SS"), ("FGHI", "")]


def test_constraint_resolver_valid_constraint():
    seq, constr = constraint_resolver("ACDE", "SS")
    assert seq == "ACDE"
    assert constr == "SS"


def test_constraint_resolver_partial_constraint_found():
    seq, constr = constraint_resolver("ACDE", "SCNT")
    assert seq == "ACDE"
    assert constr == "SCNT"


def test_constraint_resolver_partial_constraint_fallback():
    seq, constr = constraint_resolver("ACDE", "SC")
    assert seq == "ACDE"
    assert constr in {"SCNT", "SCCT", "SCSC"}  # any of the SC* matches


def test_constraint_resolver_empty_constraint():
    seq, constr = constraint_resolver("ACDE", "")
    assert seq == "ACDE"
    assert constr == ""


def test_constraint_resolver_invalid_constraint_raises():
    with pytest.raises(InvalidConstraintError):
        constraint_resolver("ACDE", "INVALID")


def test_process_constraints_yields_expected(tmp_path):
    fasta_content = ">seq1|SS\nACDE\n>seq2|HT\nFGHI"
    fasta_file = tmp_path / "test2.fasta"
    fasta_file.write_text(fasta_content)

    results = list(process_constraints(fasta_file))
    assert results == [("ACDE", "SS"), ("FGHI", "HT")]
