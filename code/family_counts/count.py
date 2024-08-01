import csv
from collections import Counter
import sys
import os

def count_mirna_families(input_file):
    # Read the TSV file and count miRNA families
    mirna_families = []
    with open(input_file, 'r') as file:
        tsv_reader = csv.DictReader(file, delimiter='\t')
        for row in tsv_reader:
            mirna_families.append(row['noncodingRNA_fam'])

    # Count occurrences of each miRNA family
    family_counts = Counter(mirna_families)

    # Sort the results from highest to lowest frequency
    sorted_families = sorted(family_counts.items(), key=lambda x: x[1], reverse=True)

    return sorted_families

def write_results_to_tsv(sorted_families, output_file):
    with open(output_file, 'w', newline='') as file:
        tsv_writer = csv.writer(file, delimiter='\t')
        tsv_writer.writerow(['miRNA Family', 'Occurrences'])
        tsv_writer.writerows(sorted_families)

def main():
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <path_to_input_tsv_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    
    try:
        sorted_families = count_mirna_families(input_file)

        # Generate output file name
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file = f"{base_name}_mirna_family_counts.tsv"

        write_results_to_tsv(sorted_families, output_file)

        print(f"Results have been written to {output_file}")
        print(f"Total unique miRNA families: {len(sorted_families)}")
    
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)
    except KeyError:
        print("Error: The TSV file does not contain a 'noncodingRNA_fam' column.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()