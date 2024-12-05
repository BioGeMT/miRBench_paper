import argparse
import csv
from collections import defaultdict
import os

def count_column_values(filename, count_column):
    value_counts = defaultdict(int)
    with open(filename, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            value = row[count_column]
            value_counts[value] += 1
    return value_counts

def parse_arguments():
    parser = argparse.ArgumentParser(description='Count occurrences in TSV files based on specified column.')
    parser.add_argument('input_files', nargs=3, help='Input TSV files')
    parser.add_argument('output_file', help='Output TSV file')
    parser.add_argument('count_column', help='Name of the column to count values from')
    return parser.parse_args()

def process_input_files(input_files, count_column):
    all_counts = []
    all_values = set()
    for input_file in input_files:
        counts = count_column_values(input_file, count_column)
        all_counts.append(counts)
        all_values.update(counts.keys())
    return all_counts, all_values

def get_column_names(input_files, count_column):
    return [count_column, 'Manakov2022', 'Hejret2023', 'Klimentova2022']

def sort_values(all_values, all_counts):
    return sorted(all_values, 
                 key=lambda x: sum(counts.get(x, 0) for counts in all_counts), 
                 reverse=True)

def write_output(output_file, column_names, sorted_values, all_counts):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(column_names)
        for value in sorted_values:
            row = [value] + [counts.get(value, 0) for counts in all_counts]
            writer.writerow(row)

def main():
    args = parse_arguments()
    all_counts, all_values = process_input_files(args.input_files, args.count_column)
    column_names = get_column_names(args.input_files, args.count_column)
    sorted_values = sort_values(all_values, all_counts)
    write_output(args.output_file, column_names, sorted_values, all_counts)

if __name__ == "__main__":
    main()