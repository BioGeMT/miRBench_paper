from miRBench.encoder import get_encoder
from miRBench.predictor import get_predictor
from miRBench.dataset import get_dataset_df
import pandas as pd
import argparse

def get_seeds(df):
    seed_types = ["Seed6mer", "Seed6merBulgeOrMismatch"]
    for tool in seed_types:       
        encoder = get_encoder(tool)
        predictor = get_predictor(tool)
        encoded_input = encoder(df)
        output = predictor(encoded_input)
        df[tool] = output
    return df
        
def main():
    parser = argparse.ArgumentParser(description="Add seed types via miRBench")
    parser.add_argument("--ifile", type=str, help="Input file")
    parser.add_argument("--ofile", type=str, help="Output file with seed types")
    args = parser.parse_args()

    # Read input file
    df = pd.read_csv(args.ifile, sep='\t')

    # Add seed types
    df_seedtypes = get_seeds(df)

    # Write seed types to file
    df_seedtypes.to_csv(args.ofile, sep='\t', index=False)

if __name__ == "__main__":
    main()