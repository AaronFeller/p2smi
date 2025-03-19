from p2smi.genPeps import (
    CONSTRAINTS,
    build_sequence,
    calculate_amino_acid_counts,
    generate_sequences,
    get_amino_acid_lists,
)


def test_calculate_amino_acid_counts_sum():
    seq_len = 50
    ln_c, dn_c, ln_nc, dn_nc = calculate_amino_acid_counts(seq_len, 0.2, 0.3)
    assert ln_c + dn_c + ln_nc + dn_nc == seq_len


def test_build_sequence_length_and_content():
    amino_lists = get_amino_acid_lists()
    seq = build_sequence(30, 0.2, 0.2, amino_lists)
    assert len(seq) == 30


def test_generate_sequences_count_and_constraints():
    sequences = generate_sequences(5, 10, 12, 0.2, 0.2, constraints=["SS"])
    assert len(sequences) == 5
    for seq_id, seq in sequences.items():
        assert seq_id.startswith("seq_")
        assert "|SS" in seq_id
        assert 10 <= len(seq) <= 12


def test_generate_sequences_zero_length():
    sequences = generate_sequences(2, 0, 0, 0.5, 0.5, constraints=[])
    for seq in sequences.values():
        assert seq == ""


def test_CONSTRAINTS():
    assert isinstance(CONSTRAINTS, list)
    assert all(
        isinstance(constraint, str) for constraint in CONSTRAINTS
    )  # Check if all constraints are strings
    assert len(CONSTRAINTS) > 0  # Ensure there are some constraints defined
