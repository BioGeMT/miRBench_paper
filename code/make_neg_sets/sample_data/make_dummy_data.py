import random
import string
import pandas as pd

# Function to generate a random DNA sequence of a given length
def generate_random_sequence(length):
    return ''.join(random.choices('ACGT', k=length))

# Generate dummy data
def generate_dummy_data(num_samples):
    data = []
    noncodingRNA_fams = [f'mir-{i}' for i in range(1, 301)]
    features = ['three_prime_utr', 'five_prime_utr', 'intron', 'exon']
    
    # Create a mapping of noncodingRNA_fams to unique seq_m
    seq_m_mapping = {fam: generate_random_sequence(22) for fam in noncodingRNA_fams}
    
    for i in range(num_samples):
        seq_g = generate_random_sequence(50)
        noncodingRNA_fam = noncodingRNA_fams[i % len(noncodingRNA_fams)]
        seq_m = seq_m_mapping[noncodingRNA_fam]
        feature = random.choice(features)
        test = random.choice([True, False])
        label = 1  # Since we are generating positive samples
        data.append([seq_g, seq_m, noncodingRNA_fam, feature, test, label])
    
    return data

# Number of samples to generate
num_samples = 1000

# Generate the data
dummy_data = generate_dummy_data(num_samples)

# Create a DataFrame and save to a file
df = pd.DataFrame(dummy_data, columns=['seq.g', 'seq.m', 'noncodingRNA_fam', 'feature', 'test', 'label'])
df.to_csv('dummy_data.tab', sep='\t', index=False)

print("Dummy data generated and saved to 'dummy_data.tab'")
