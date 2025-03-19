import argparse
import random

from p2smi.utilities.aminoacids import all_aminos


def get_amino_acid_lists():
    all_aas = [letter for entry in all_aminos.values() for letter in entry["Letter"]]
    canonical = [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    ]
    lower_canon = [aa.lower() for aa in canonical if aa != "G"]
    upper_noncanonical = [
        aa
        for aa in all_aas
        if aa not in canonical and aa not in lower_canon and aa.isupper()
    ]
    lower_noncanonical = [
        aa
        for aa in all_aas
        if aa not in canonical and aa not in lower_canon and aa.islower()
    ]
    return canonical, lower_canon, upper_noncanonical, lower_noncanonical


CONSTRAINTS = ["SS", "HT", "SCNT", "SCCT", "SCSC"]


def calculate_amino_acid_counts(seq_len, noncanonical_percent, dextro_percent):
    dn_nc = round(seq_len * dextro_percent * noncanonical_percent)
    dn_c = round(seq_len * dextro_percent * (1 - noncanonical_percent))
    ln_nc = round(seq_len * (1 - dextro_percent) * noncanonical_percent)
    ln_c = seq_len - (dn_nc + dn_c + ln_nc)
    return ln_c, dn_c, ln_nc, dn_nc


def build_sequence(seq_len, noncanonical_percent, dextro_percent, amino_lists):
    canonical, lower_canon, upper_noncon, lower_noncon = amino_lists
    ln_c, dn_c, ln_nc, dn_nc = calculate_amino_acid_counts(
        seq_len, noncanonical_percent, dextro_percent
    )
    sequence_parts = (
        random.choices(canonical, k=ln_c)
        + random.choices(lower_canon, k=dn_c)
        + random.choices(upper_noncon, k=ln_nc)
        + random.choices(lower_noncon, k=dn_nc)
    )
    random.shuffle(sequence_parts)
    return "".join(sequence_parts)


def generate_sequences(
    num_sequences,
    min_length,
    max_length,
    noncanonical_percent,
    dextro_percent,
    constraints,
):
    amino_lists = get_amino_acid_lists()

    def make_sequence(i):
        seq_id = f"seq_{i + 1}"
        constraint = random.choice(constraints) if constraints else None
        if constraint:
            seq_id += f"|{constraint}"
        seq_len = random.randint(min_length, max_length)
        seq = build_sequence(seq_len, noncanonical_percent, dextro_percent, amino_lists)
        return seq_id, seq

    return dict(make_sequence(i) for i in range(num_sequences))


def output_sequences(sequences, outfile=None):
    lines = [f">{seq_id}\n{seq}" for seq_id, seq in sequences.items()]
    output = "\n".join(lines)
    if outfile:
        with open(outfile, "w") as f:
            f.write(output + "\n")
    else:
        print(output)


def main():
    parser = argparse.ArgumentParser(
        description="Generate random amino acid sequences."
    )
    parser.add_argument("--min_seq_len", type=int, default=5)
    parser.add_argument("--max_seq_len", type=int, default=100)
    parser.add_argument("--noncanonical_percent", type=float, default=0.2)
    parser.add_argument("--lowercase_percent", type=float, default=0.2)
    parser.add_argument("--num_sequences", type=int, default=1)
    parser.add_argument("--constraints", type=str, nargs="*", default=[])
    parser.add_argument("--outfile", type=str, default=None)
    args = parser.parse_args()

    constraints = (
        CONSTRAINTS
        if args.constraints == ["all"] or not args.constraints
        else args.constraints
    )

    sequences = generate_sequences(
        args.num_sequences,
        args.min_seq_len,
        args.max_seq_len,
        args.noncanonical_percent,
        args.lowercase_percent,
        constraints,
    )

    output_sequences(sequences, args.outfile)


if __name__ == "__main__":
    main()
