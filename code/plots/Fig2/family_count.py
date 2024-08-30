import argparse
import csv
from collections import defaultdict
import os

def count_mirna_families(filename):
    family_counts = defaultdict(int)
    with open(filename, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            family = row['miRNA_fam']
            family_counts[family] += 1
    return family_counts

def parse_arguments():
    parser = argparse.ArgumentParser(description='Count miRNA families in TSV files.')
    parser.add_argument('input_files', nargs=3, help='Input TSV files')
    parser.add_argument('output_file', help='Output TSV file')
    return parser.parse_args()

def process_input_files(input_files):
    all_counts = []
    all_families = set()
    for input_file in input_files:
        counts = count_mirna_families(input_file)
        all_counts.append(counts)
        all_families.update(counts.keys())
    return all_counts, all_families

def get_column_names(input_files):
    return ['miRNA Family'] + [os.path.splitext(os.path.basename(f))[0] for f in input_files]

def sort_families(all_families, all_counts):
    return sorted(all_families, key=lambda x: all_counts[0].get(x, 0), reverse=True)

def write_output(output_file, column_names, sorted_families, all_counts):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(column_names)
        for family in sorted_families:
            row = [family] + [counts.get(family, 0) for counts in all_counts]
            writer.writerow(row)
    print(f"Results written to {output_file}")

def main():
    args = parse_arguments()
    all_counts, all_families = process_input_files(args.input_files)
    column_names = get_column_names(args.input_files)
    sorted_families = sort_families(all_families, all_counts)
    write_output(args.output_file, column_names, sorted_families, all_counts)

if __name__ == "__main__":
    main()