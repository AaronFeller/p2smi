# CycloPs_v2: A Cyclic Peptide Library Generator

**CycloPs_v2** is a cyclic peptide library generator designed to facilitate the creation and exploration of diverse peptide libraries, including noncanonical and cyclic peptides. This package builds on interest generated by the **PeptideCLM model**, as described in the [bioRxiv paper](https://www.biorxiv.org/content/10.1101/2024.08.09.607221v1).

> **Please cite the PeptideCLM paper** (or its published version) if you use this package. A JOSS paper on CycloPs_v2 is forthcoming.

## Usage Guide

1. **Generate Peptides**:
   - Run `generate_randpeps.py` to generate random peptide sequences.
   
2. **Generate SMILES Notation**:
   - Pass the resulting FASTA file to `PepLibGen/StructGen/fasta.py` for SMILES generation. 
   - The output is a text file, with each line containing:
     ```
     {Comma-separated Amino Acids}-{cyclization notation}: {SMILES string}
     ```

3. **Noncanonical Modifications**:
   - The Jupyter notebook `build_noncanon_aa.ipynb` provides workflows for adding modifications like N-methylation and PEGylation.
   - These modifications will be packaged as functions in future releases.

4. **Amino Acid Notation**:
   - Notation for single-character amino acids is found in `PepLibGen/StructGen/aminoacids.py`.

## Future Work
Improvements, including modular functions for peptide modifications and additional SMILES generation options, are in progress.

---

CycloPs_v2 is a work in progress. Feedback and contributions are welcome as we continue to develop and refine this tool.
