import pandas as pd
import argparse
import os
from Bio.Seq import Seq

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_seed_match(target, mirna):
    rc_target = reverse_complement(target)
    seed6 = mirna[1:7]
    seed7 = mirna[1:8]
    seed8 = mirna[1:9]
    
    if seed8 in rc_target:
        return '8mer'
    elif seed7 in rc_target:
        return '7mer'
    elif seed6 in rc_target:
        return '6mer'
    else:
        return 'none'

def main(input_file):
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = f"{base_name}_verification_report.txt"

    df = pd.read_csv(input_file, sep='\t')

    with open(output_file, 'w') as report:
        report.write("Verification Report for Seed Match Functions\n\n")

        report.write("1. Testing reverse_complement function:\n")
        for i in range(10):  # Test with 10 entries
            test_seq = df['seq.g'].iloc[i]
            rc_seq = reverse_complement(test_seq)
            report.write(f"Test {i+1}:\n")
            report.write(f"  Input sequence:     {test_seq}\n")
            report.write(f"  Reverse complement: {rc_seq}\n")
            report.write(f"  Double RC:          {reverse_complement(rc_seq)}\n")
            report.write(f"  Test passed: {test_seq == reverse_complement(rc_seq)}\n\n")

        report.write("2. Testing find_seed_match function:\n")
        for i in range(10):  # Test with 10 entries
            target = df['seq.g'].iloc[i]
            mirna = df['seq.m'].iloc[i]
            result = find_seed_match(target, mirna)
            rc_target = reverse_complement(target)
            report.write(f"Test {i+1}:\n")
            report.write(f"  Target:    {target}\n")
            report.write(f"  miRNA:     {mirna}\n")
            report.write(f"  RC Target: {rc_target}\n")
            report.write(f"  6mer seed: {mirna[1:7]}\n")
            report.write(f"  7mer seed: {mirna[1:8]}\n")
            report.write(f"  8mer seed: {mirna[1:9]}\n")
            report.write(f"  Result:    {result}\n")
            report.write(f"  6mer in RC: {'Yes' if mirna[1:7] in rc_target else 'No'}\n")
            report.write(f"  7mer in RC: {'Yes' if mirna[1:8] in rc_target else 'No'}\n")
            report.write(f"  8mer in RC: {'Yes' if mirna[1:9] in rc_target else 'No'}\n\n")

        report.write("3. Analyzing seed match occurrences in the entire dataset:\n")
        counts = {'8mer': 0, '7mer': 0, '6mer': 0, 'none': 0}
        for _, row in df.iterrows():
            match_type = find_seed_match(row['seq.g'], row['seq.m'])
            counts[match_type] += 1
        total = sum(counts.values())
        for match_type, count in counts.items():
            percentage = (count / total) * 100
            report.write(f"{match_type}: {count} ({percentage:.2f}%)\n")
        report.write("\n")

        report.write("4. Verifying consistency of match types:\n")
        inconsistencies = 0
        for _, row in df.iterrows():
            target = row['seq.g']
            mirna = row['seq.m']
            rc_target = reverse_complement(target)
            match_type = find_seed_match(target, mirna)
            if match_type == '8mer' and mirna[1:9] not in rc_target:
                inconsistencies += 1
            elif match_type == '7mer' and mirna[1:8] not in rc_target:
                inconsistencies += 1
            elif match_type == '6mer' and mirna[1:7] not in rc_target:
                inconsistencies += 1
            elif match_type == 'none' and (mirna[1:7] in rc_target or mirna[1:8] in rc_target or mirna[1:9] in rc_target):
                inconsistencies += 1
        report.write(f"Total inconsistencies found: {inconsistencies}\n")
        report.write(f"Percentage of inconsistencies: {(inconsistencies/len(df))*100:.2f}%\n")

    print(f"Verification complete. Results written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verification of seed matching functions")
    parser.add_argument("input_file", help="Path to the Klimentova TSV file")
    args = parser.parse_args()
    main(args.input_file)

