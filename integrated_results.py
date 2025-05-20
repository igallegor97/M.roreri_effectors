import os
import glob
import csv
from Bio import SeqIO
import matplotlib.pyplot as plt
import re

def load_lengths(fasta_files):
    lengths = {}
    for file in fasta_files:
        for record in SeqIO.parse(file, "fasta"):
            lengths[record.id] = len(record.seq)
    return lengths

def load_signalp_predictions(gff_files):
    signalp_hits = set()
    for file in gff_files:
        print(f"\nLeyendo archivo SignalP: {file}")  # Debug
        with open(file) as f:
            for line in f:
                if not line.startswith("#") and line.strip():
                    cols = line.strip().split('\t')
                    print(f"Línea procesada: {cols}")  # Debug
                    if len(cols) >= 4 and cols[2] == "signal_peptide":  # Cambié cols[3] → cols[2] (índice base 0)
                        protein_id = cols[0]
                        print(f"SignalP hit encontrado: {protein_id}")  # Debug
                        signalp_hits.add(protein_id)
    print(f"\nTotal SignalP hits: {len(signalp_hits)}")  # Debug
    return signalp_hits

def load_tmhmm_predictions(gff_files):
    tmhmm_hits = set()
    for file in gff_files:
        print(f"\nLeyendo archivo TMHMM: {file}")  # Debug
        with open(file) as f:
            for line in f:
                if line.startswith("#") and "Number of predicted TMRs:" in line:
                    parts = line.strip().split()
                    protein_id = parts[1]  # Asumiendo formato "# tr_ID_ID_MONRR"
                    num_tmrs = int(parts[-1])
                    print(f"Proteína: {protein_id}, TMRs: {num_tmrs}")  # Debug
                    if num_tmrs > 0:
                        tmhmm_hits.add(protein_id)
    print(f"\nTotal TMHMM hits: {len(tmhmm_hits)}")  # Debug
    return tmhmm_hits

def load_effector_predictions(fasta_files):
    effector_ids = set()
    for file in fasta_files:
        for record in SeqIO.parse(file, "fasta"):
            effector_ids.add(record.id)
    return effector_ids

def plot_summary(total, n_signal, n_tm, n_effector, n_final):
    labels = [
        'Total proteins',
        'Signal peptide',
        'Transmembrane domains',
        'EffectorP effectors',
        'Final candidates'
    ]
    values = [total, n_signal, n_tm, n_effector, n_final]
    colors = ['#d9d9d9', '#80b1d3', '#fb8072', '#b3de69', '#ffed6f']

    plt.figure(figsize=(10, 6))
    bars = plt.bar(labels, values, color=colors)
    plt.title('Summary of Protein Features')
    plt.ylabel('Number of proteins')
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()

    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2.0, yval + 10, int(yval), ha='center', va='bottom')

    plt.savefig("summary_plot.png")
    print("Summary plot saved as summary_plot.png")

def integrate_data(lengths, signalp_hits, tmhmm_hits, effector_ids, output_csv, filtered_csv, summary_file=None):
    total = len(lengths)
    n_signal = 0
    n_tm = 0
    n_effector = 0
    n_final_candidates = 0

    with open(output_csv, "w", newline='') as out_all, open(filtered_csv, "w", newline='') as out_filtered:
        writer_all = csv.writer(out_all)
        writer_filtered = csv.writer(out_filtered)

        header = ["Protein ID", "Signal Peptide", "Transmembrane Domain", "EffectorP Prediction", "Length (aa)"]
        writer_all.writerow(header)
        writer_filtered.writerow(header)

        for protein_id in lengths:
            signal = "Yes" if protein_id in signalp_hits else "No"
            tm = "Yes" if protein_id in tmhmm_hits else "No"
            effector = "Effector" if protein_id in effector_ids else "Non-effector"
            length = lengths[protein_id]

            if signal == "Yes":
                n_signal += 1
            if tm == "Yes":
                n_tm += 1
            if effector == "Effector":
                n_effector += 1

            row = [protein_id, signal, tm, effector, length]
            writer_all.writerow(row)

            if signal == "Yes" and tm == "No" and effector == "Effector":
                writer_filtered.writerow(row)
                n_final_candidates += 1

    summary = f"""
Summary of Results
==================
Total proteins analyzed:          {total}
Proteins with signal peptide:     {n_signal}
Proteins with TM domains:         {n_tm}
Proteins predicted as effectors:  {n_effector}
Final secreted effector candidates (SignalP+EffectorP - TMHMM): {n_final_candidates}
"""
    print(summary)
    if summary_file:
        with open(summary_file, "w") as f:
            f.write(summary)

    plot_summary(total, n_signal, n_tm, n_effector, n_final_candidates)

if __name__ == "__main__":
    base_dir = "."

    fasta_inputs = glob.glob(os.path.join(base_dir, "moniliophthora_filtered_split_*.fasta"))
    signalp_gffs = glob.glob(os.path.join(base_dir, "signal_results*", "*.gff3"))
    tmhmm_gffs = glob.glob(os.path.join(base_dir, "TMRs_*.gff3"))
    effector_fastas = glob.glob(os.path.join(base_dir, "EffectorCandidates_*.fasta"))

    print("Reading sequences and lengths...")
    lengths = load_lengths(fasta_inputs)

    print("Reading SignalP predictions...")
    signalp_hits = load_signalp_predictions(signalp_gffs)

    print("Reading TMHMM predictions...")
    tmhmm_hits = load_tmhmm_predictions(tmhmm_gffs)

    print("Reading EffectorP predictions...")
    effector_ids = load_effector_predictions(effector_fastas)

    print("Integrating results and generating plot...")
    integrate_data(
        lengths,
        signalp_hits,
        tmhmm_hits,
        effector_ids,
        output_csv="integrated_results.csv",
        filtered_csv="filtered_effectors.csv",
        summary_file="summary.txt"
    )

    print("Done")
