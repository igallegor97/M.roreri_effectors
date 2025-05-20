# Integration Script: Prediction of Secreted Effectors in *Moniliophthora roreri*

This script, `integrate_results.py`, performs the integration and filtering of prediction results obtained from three bioinformatics tools:

- **SignalP 6.0**: predicts signal peptides (secretion signals).
- **DeepTMHMM**: predicts transmembrane domains.
- **EffectorP 3.0**: predicts whether a protein is a potential effector.

The main goal is to identify secreted effector proteins without transmembrane domains, using outputs from the above tools and a FASTA file containing protein sequences.

---

## What Does This Script Do?

1. **Reads protein sequences** and calculates their lengths.
2. **Parses SignalP outputs** (GFF3 format) to find proteins with signal peptides.
3. **Parses DeepTMHMM outputs** (GFF3 format) to detect proteins with TM domains.
4. **Parses EffectorP outputs** (FASTA format) to identify predicted effectors.
5. **Combines all results into a table** and applies filtering:
    - Keeps proteins with signal peptides (`SignalP = Yes`),
    - Without TM domains (`TMHMM = No`),
    - And predicted as effectors (`EffectorP = Effector`).
6. **Generates a summary report** and a **bar plot**.

## Required Input Files

These should be generated beforehand using the respective tools:

| Tool        | Input Type           | Example Pattern                    |
|-------------|----------------------|------------------------------------|
| SignalP     | `.gff3`              | `signal_results*/results.gff3`     |
| DeepTMHMM   | `.gff3`              | `TMRs_*.gff3`                       |
| EffectorP   | `.fasta`             | `EffectorCandidates_*.fasta`       |
| Sequences   | `.fasta` (original)  | `moniliophthora_filtered_split_*.fasta` |

## Output Files

| File Name                  | Description                                           |
|---------------------------|-------------------------------------------------------|
| `integrated_results.csv`  | All proteins with their predictions.                  |
| `filtered_effectors.csv`  | Only those matching final filter (signal + effector - TM). |
| `summary.txt`             | Text summary of counts.                               |
| `summary_plot.png`        | Bar chart of overall results.                         |

 License and Credits
This script was developed as part of the project:

> Prediction of Secreted Effectors in Moniliophthora roreri
> Bioinformatica e Genomica, cen5789
> CENA, Universidade de São Paulo, 2025

Developed by Isabella Gallego Rendón.
