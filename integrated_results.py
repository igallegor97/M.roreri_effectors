import os
import glob
import csv
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib_venn import venn3

def load_lengths(fasta_files):
    """Loads protein lengths from FASTA files."""
    lengths = {}
    for file in fasta_files:
        for record in SeqIO.parse(file, "fasta"):
            lengths[record.id] = len(record.seq)
    return lengths

def load_signalp_predictions(gff_files):
    """Extracts proteins with signal peptides from SignalP GFF files."""
    signalp_hits = set()
    for file in gff_files:
        with open(file) as f:
            for line in f:
                if not line.startswith("#") and line.strip():
                    cols = line.strip().split('\t')
                    if len(cols) >= 4 and cols[2] == "signal_peptide":
                        signalp_hits.add(cols[0])
    return signalp_hits

def load_tmhmm_predictions(gff_files):
    """Extracts proteins with transmembrane domains from TMHMM GFF files."""
    tmhmm_hits = set()
    for file in gff_files:
        with open(file) as f:
            for line in f:
                if line.startswith("#") and "Number of predicted TMRs:" in line:
                    parts = line.strip().split()
                    if int(parts[-1]) > 0:
                        tmhmm_hits.add(parts[1])
    return tmhmm_hits

def load_effector_predictions(fasta_files):
    """Loads IDs of proteins predicted as effectors by EffectorP."""
    return {record.id for file in fasta_files for record in SeqIO.parse(file, "fasta")}

def plot_summary(total, n_signal, n_tm, n_effector, n_final):
    """Summary bar chart."""
    labels = ['Total proteins', 'Signal peptide', 'TM domains', 'EffectorP effectors', 'Final candidates']
    values = [total, n_signal, n_tm, n_effector, n_final]
    colors = ['#d9d9d9', '#80b1d3', '#fb8072', '#b3de69', '#ffed6f']

    plt.figure(figsize=(10, 6))
    bars = plt.bar(labels, values, color=colors)
    plt.title('Summary of Protein Features')
    plt.ylabel('Count')
    plt.xticks(rotation=30, ha='right')
    
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 5, int(yval), ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig("summary_plot.png")

def plot_length_distribution(lengths, final_candidates):
    """Distribution of lengths: all proteins vs final candidates."""
    plt.figure(figsize=(10, 6))
    plt.hist([list(lengths.values()), [lengths[pid] for pid in final_candidates]],
             bins=30, alpha=0.7, label=['All proteins', 'Final candidates'])
    plt.xlabel('Protein Length (aa)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.title('Protein Length Distribution')
    plt.savefig("length_distribution.png")

def plot_venn_diagram(signalp, tmhmm, effectors):
    """Venn diagram showing overlap between tools."""
    plt.figure(figsize=(8, 8))
    venn3([signalp, effectors, tmhmm], ('SignalP', 'EffectorP', 'TMHMM'))
    plt.title("Tool Prediction Overlap")
    plt.savefig("venn_diagram.png")

def plot_feature_correlation(lengths, signalp, tmhmm, effectors):
    """Heatmap showing correlation between features."""
    data = []
    for pid in lengths:
        data.append({
            'Length': lengths[pid],
            'SignalP': int(pid in signalp),
            'TMHMM': int(pid in tmhmm),
            'EffectorP': int(pid in effectors)
        })
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(pd.DataFrame(data).corr(), annot=True, cmap='coolwarm', vmin=-1, vmax=1)
    plt.title("Feature Correlation")
    plt.tight_layout()
    plt.savefig("correlation_heatmap.png")

def integrate_data(lengths, signalp, tmhmm, effectors, output_csv, filtered_csv, summary_file):
    """Integrates data and generates plots."""
    total = len(lengths)
    stats = {
        'signal': sum(1 for pid in lengths if pid in signalp),
        'tm': sum(1 for pid in lengths if pid in tmhmm),
        'effector': sum(1 for pid in lengths if pid in effectors),
        'final': []
    }

    with open(output_csv, 'w') as f_all, open(filtered_csv, 'w') as f_filtered:
        writer_all = csv.writer(f_all)
        writer_filtered = csv.writer(f_filtered)
        header = ["Protein ID", "SignalP", "TMHMM", "EffectorP", "Length"]
        writer_all.writerow(header)
        writer_filtered.writerow(header)

        for pid in lengths:
            row = [pid, 
                   'Yes' if pid in signalp else 'No',
                   'Yes' if pid in tmhmm else 'No',
                   'Effector' if pid in effectors else 'Non-effector',
                   lengths[pid]]
            writer_all.writerow(row)
            
            if pid in signalp and pid not in tmhmm and pid in effectors:
                writer_filtered.writerow(row)
                stats['final'].append(pid)

    # Plot generation
    plot_summary(total, stats['signal'], stats['tm'], stats['effector'], len(stats['final']))
    plot_length_distribution(lengths, stats['final'])
    plot_venn_diagram(signalp, tmhmm, effectors)
    plot_feature_correlation(lengths, signalp, tmhmm, effectors)

    # Text summary
    summary = f"""Summary of Results:
Total proteins: {total}
With signal peptide: {stats['signal']}
With TM domains: {stats['tm']}
EffectorP effectors: {stats['effector']}
Final candidates: {len(stats['final'])}"""
    
    with open(summary_file, 'w') as f:
        f.write(summary)
    print(summary)

if __name__ == "__main__":
    # Path configuration
    base_dir = "."
    files = {
        'fasta': glob.glob(os.path.join(base_dir, "moniliophthora_filtered_split_*.fasta")),
        'signalp': glob.glob(os.path.join(base_dir, "signal_results*", "*.gff3")),
        'tmhmm': glob.glob(os.path.join(base_dir, "TMRs_*.gff3")),
        'effector': glob.glob(os.path.join(base_dir, "EffectorCandidates_*.fasta"))
    }

    # Data loading
    print("Loading data")
    data = {
        'lengths': load_lengths(files['fasta']),
        'signalp': load_signalp_predictions(files['signalp']),
        'tmhmm': load_tmhmm_predictions(files['tmhmm']),
        'effectors': load_effector_predictions(files['effector'])
    }

    # Processing and integration
    print("Generating outputs...")
    integrate_data(
        data['lengths'],
        data['signalp'],
        data['tmhmm'],
        data['effectors'],
        "integrated_results.csv",
        "filtered_effectors.csv",
        "summary.txt"
    )
    print("Done.")
