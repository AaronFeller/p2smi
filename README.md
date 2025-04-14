# p2smi: Peptide FASTA-to-SMILES Conversion and Molecular Property Tools

**p2smi** is a Python package for generating and modifying peptide SMILES strings from FASTA input and computing molecular properties. It supports cyclic and linear peptides, noncanonical amino acids, and common chemical modifications (e.g., N-methylation, PEGylation).

This package was released in its current form to support work on the **PeptideCLM** model, described in our [Publication](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01441).

> **If you use this tool, please cite the PeptideCLM paper.** A publication for this specific toolkit is forthcoming.

## Manuscript
- [View PDF](manuscript/paper.pdf)
- [View markdown source](manuscript/paper.md)

## Directory

- [Features](#features)
- [Installation](#installation)
- [Command-Line Tools](#command-line-tools)
- [Example Usage](#example-usage)
- [Future Work](#future-work)
- [Contributing](#for-contributors)
- [License](#license)
- [Citation](#citation)

## Features
- Convert peptide FASTA files into valid SMILES strings
- Automatically handle peptide cyclizations (disulfide, head-to-tail, side-chain to N-term, side-chain to C-term, side-chain to side-chain)
- Modify peptide SMILES with customizable N-methylation and PEGylation
- Evaluate synthesis feasibility with defined synthesis rules
- Compute molecular properties: logP, TPSA, molecular formula, and Lipinski rule evaluation

## Installation
```bash
pip install p2smi
```
For development:
```bash
git clone https://github.com/AaronFeller/p2smi.git
cd p2smi
pip install -e .[dev]
```

## Command-Line Tools

| Command               | Description                                                     |
|-----------------------|-----------------------------------------------------------------|
| `generate-peptides`  | Generate random peptide sequences based on user-defined constraints and modifications |
| `fasta2smi`          | Convert a FASTA file of peptide sequences into SMILES format    |
| `modify-smiles`      | Apply modifications (N-methylation, PEGylation) to existing SMILES strings |
| `smiles-props`       | Compute molecular properties (logP, TPSA, formula, Lipinski rules) from SMILES |
| `synthesis-check`    | Check synthesis constraints for peptides (*currently only functional for natural amino acids*) |

> Run each command with `--help` to view usage and options:
```bash
generate-peptides --help
fasta2smi --help
modify-smiles --help
smiles-props --help
synthesis-check --help
```

## Example Usage

**Generate a random peptide:**
```bash
generate-peptides \
    --max_seq_len 20 \
    --min_seq_len 10 \
    --noncanonical_percent 0.1 \
    --lowercase_percent 0.1 \
    --num_sequences 10 \
    --constraints all \
    --outfile outputfile.smi
```

**Convert a FASTA file to SMILES:**
```bash
fasta2smi -i peptides.fasta -o output.smi
```

**Modify existing SMILES strings (N-methylation/PEGylation):**
```bash
modify-smiles -i input.smi -o modified.smi --peg_rate 0.3 --nmeth_rate 0.2 --nmeth_residues 0.25
```

**Compute properties of a SMILES string:**
```bash
smiles-props "C1CC(NC(=O)C2CC2)C1"
```

**Check synthetic feasability**
```bash
synthesis-check output.smi  # only works for natural amino acids
```

## Future Work
- Expand support for additional post-translational modifications (build importer)
- Enhance synthesis-check with rules for noncanonical amino acid and modified peptides
- Expand usage of mol files (applying RDKit's Chem.MolFromSmiles() function)
- Include alternative encodings (HELM, SELFIES, etc.)
- Enable batch processing/threading for high-throughput analysis
- Incorporate predictive models for synthesis of unnatural amino acids

## For Contributors
There are several ways you can contribute to this project:

- Reporting Bugs: If you encounter any issues or unexpected behavior, please let us know by opening an issue.
- Suggesting Enhancements: Have ideas to improve the project? We’d love to hear them! Share your suggestions by opening an issue.
- Submitting Pull Requests: If you’d like to fix a bug or implement a new feature, you can submit a pull request.
- Improving Documentation: Clear and comprehensive documentation helps everyone.

## License
[MIT License](https://github.com/AaronFeller/p2smi/blob/master/LICENSE)

## Citation
> If you use this tool, please cite:  
- [Peptide-Aware Chemical Language Model Successfully Predicts Membrane Diffusion of Cyclic Peptides (JCIM)](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01441)  
A JOSS paper will follow.
