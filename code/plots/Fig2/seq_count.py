import argparse
import csv
from collections import defaultdict
import os

def count_noncoding_sequences(filename):
    sequence_counts = defaultdict(int)
    with open(filename, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sequence = row['noncodingRNA']
            sequence_counts[sequence] += 1
    return sequence_counts

def parse_arguments():
    parser = argparse.ArgumentParser(description='Count unique noncodingRNA sequences in TSV files.')
    parser.add_argument('input_files', nargs=3, help='Input TSV files')
    parser.add_argument('output_file', help='Output TSV file')
    return parser.parse_args()

def process_input_files(input_files):
    all_counts = []
    all_sequences = set()
    for input_file in input_files:
        counts = count_noncoding_sequences(input_file)
        all_counts.append(counts)
        all_sequences.update(counts.keys())
    return all_counts, all_sequences

def get_column_names(input_files):
    # Use fixed column names as shown in the example
    return ['noncodingRNA_sequence', 'Manakov2022', 'Hejret2023', 'Klimentova2022']

def sort_sequences(all_sequences, all_counts):
    # Sort by total count across all files (descending)
    return sorted(all_sequences, 
                 key=lambda x: sum(counts.get(x, 0) for counts in all_counts), 
                 reverse=True)

def write_output(output_file, column_names, sorted_sequences, all_counts):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(column_names)
        for sequence in sorted_sequences:
            row = [sequence] + [counts.get(sequence, 0) for counts in all_counts]
            writer.writerow(row)
    print(f"Results written to {output_file}")

def main():
    args = parse_arguments()
    all_counts, all_sequences = process_input_files(args.input_files)
    column_names = get_column_names(args.input_files)
    sorted_sequences = sort_sequences(all_sequences, all_counts)
    write_output(args.output_file, column_names, sorted_sequences, all_counts)

if __name__ == "__main__":
    main()
