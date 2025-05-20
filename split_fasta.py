from Bio import SeqIO

input_file = "moniliophthora_proteome.fasta"
output_prefix = "moniliophthora_filtered_split"
max_seqs = 250
min_length = 10

# Read and filter sequences 

records = [record for record in SeqIO.parse(input_file, "fasta") if len(record.seq) >= min_length]

print(f"Total sequences after filtering: {len(records)}")

# Batch process 

for i in range(0, len(records), max_seqs):
    batch = records[i:i+max_seqs]
    output_file = f"{output_prefix}_{i//max_seqs + 1}.fasta"
    SeqIO.write(batch, output_file, "fasta")
    print(f"File generated: {output_file} with {len(batch)} sequences")
