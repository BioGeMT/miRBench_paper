import pandas as pd

# Load the TSV file
file_path = "Klimentova2022.tsv"  # Replace with the actual file path
data = pd.read_csv(file_path, sep='\t')

# Define the output FASTA file
fasta_path = "klimentova_sequences.fasta"

# Write the sequences from the 'gene' column to the FASTA file
with open(fasta_path, 'w') as fasta_file:
    for index, row in data.iterrows():
        sequence = row['gene']  # Access the 'gene' column
        fasta_file.write(f">Seq_{index + 1}\n{sequence}\n")

print(f"FASTA file created: {fasta_path}")
