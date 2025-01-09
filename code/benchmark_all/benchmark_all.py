from miRBench.encoder import get_encoder
from miRBench.predictor import get_predictor, list_predictors
from miRBench.dataset import get_dataset_path, list_datasets
import pandas as pd
import argparse
import os


def benchmark_all(df, dset):
    for tool in list_predictors():
        print(f"Running {tool} on {dset} dataset")           
        encoder = get_encoder(tool)
        predictor = get_predictor(tool)
        input = encoder(df)
        output = predictor(input)
        df[tool] = output
    return df


def main():
    parser = argparse.ArgumentParser(description="Benchmark all available predictors on all available datasets")
    parser.add_argument("--out_dir", type=str, default=".", help="Output directory for predictions")
    
    args = parser.parse_args()
    
    split = "test"

    # loop over all available datasets
    for dset in list_datasets():
        print(f"Downloading {dset} dataset, {split} split")
        input_file = get_dataset_path(dset, split=split)
        output_file = os.path.join(args.out_dir, f"{dset}_{split}_predictions.tsv")
        header_written = False
        df = pd.read_csv(input_file, sep='\t')
        df_preds = benchmark_all(df, dset)
        if not header_written:
            df_preds.to_csv(output_file, sep='\t', index=False, mode='w')
            header_written = True
        else:
            df_preds.to_csv(output_file, sep='\t', index=False, mode='a', header=False)

        print(f"Predictions for {dset} dataset, {splt} split, written to {output_file}")

    print("Predictions for all datasets written to the output directory")

if __name__ == "__main__":
    main()