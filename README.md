# Prediction of Secreted Effectors in *Moniliophthora roreri*

**Project for the Genomics and Bioinformatics course (cen5789, 2025)**  
**Master’s in Bioinformatics – CENA, Universidade de São Paulo**

This project focuses on the identification of secreted effector proteins in the fungal phytopathogen *Moniliophthora roreri*, the causal agent of Frosty Pod Rot in cacao.

The workflow involves:

- Splitting the proteome into manageable fragments
- Predicting signal peptides (SignalP)
- Predicting transmembrane domains (DeepTMHMM or TMHMM)
- Predicting effectors (EffectorP)
- Integrating all predictions into a final summary table

## Project Structure

```bash
moniliophthora_effectors_project/
├── data/
│ ├── moniliophthora_proteome.fasta
│ ├── moniliophthora_filtered_split_.fasta
│ ├── EffectorCandidates_.fasta
│ ├── TMRs_.gff3
│ ├── signal_results/results.gff3
├── scripts/
│ ├── split_and_filter_fasta.py
│ ├── integrate_results.py
├── results/
│ ├── integrated_results.csv
│ ├── filtered_effectors.csv
│ ├── summary.txt
│ ├── summary_plot.png
├── README.md
└── environment.yml
```

## Tools Used

| Tool          | Purpose                            | Source / Link |
|---------------|------------------------------------|---------------|
| SignalP 6.0   | Signal peptide prediction           | [DTU HealthTech](https://services.healthtech.dtu.dk/service.php?SignalP-6.0) |
| DeepTMHMM     | Transmembrane domain prediction     | [DTU HealthTech](https://services.healthtech.dtu.dk/service.php?DeepTMHMM-1.0) |
| TMHMM 2.0     | (Alternative to DeepTMHMM)          | [DTU HealthTech](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0) |
| EffectorP 3.0 | Effector prediction (fungal)        | [GitHub](https://github.com/JanaSperschneider/EffectorP-3.0) |
| Python        | Data integration, analysis          | Python ≥3.8 |
| Biopython     | Sequence parsing                    | `pip install biopython` |
| matplotlib    | Summary plotting                    | `pip install matplotlib` |

## Installation

### 1. Clone this repository

```bash
git clone https://github.com/isa-gallego/moniliophthora-effectors.git
cd moniliophthora-effectors

# Create conda env

conda create -n effector_env python=3.8
conda activate effector_env
pip install biopython matplotlib

```

## How to Run the Workflow
### 1. Split and filter the proteome
Use the provided script to split the FASTA into files of ≤250 sequences, filtering by minimum length:

```bash
python scripts/split_and_filter_fasta.py \
    --input data/moniliophthora_proteome.fasta \
    --output data/ --max_seqs 250 --min_len 10
```

### 2. Run external prediction tools
- SignalP 6.0 (CLI)

- DeepTMHMM (web)

- EffectorP 3.0 (CLI)

### 3. Integrate results

This will generate:

- integrated_results.csv: full table of all proteins

- filtered_effectors.csv: final candidate list

- summary.txt: counts of proteins by category

- summary_plot.png: bar plot of results

### Criteria for Effector Candidates
To be considered a secreted effector, a protein must:

- Have a signal peptide (SignalP = Yes)
- 

- Not contain transmembrane domains (DeepTMHMM = No)

- Be predicted as an effector (EffectorP = Effector)
