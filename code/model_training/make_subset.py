import pandas as pd
import argparse


def subsample_dataset(dataset, size, output):
    
    df = pd.read_csv(dataset, sep='\t')
    pos = df[df['label']==1]
    neg = df[df['label']==0]
    pos = pos.sample(n=size, ignore_index=True, random_state=42)
    neg = neg.sample(n=size, ignore_index=True, random_state=42)
    subsampled = pd.concat([pos, neg]).sample(frac=1, ignore_index=True, random_state=42)
    subsampled.to_csv(output, sep='\t', index=False)
    

def main():
    parser = argparse.ArgumentParser(description="Make subset of specified size out of the input tsv file")
    parser.add_argument('--N', type=int, required=True, help="Number of positive samples to be randomly chosen from the dataset. Negatives will be subsampled equally.")
    parser.add_argument('--dataset', type=str, required=True, help="File with the input tsv dataset")
    parser.add_argument('--output', type=str, required=False, help="Filename to save the subsampled dataset")
    args = parser.parse_args()

    subsample_dataset(args.dataset, args.N, args.output)

if __name__ == "__main__":
    main()
