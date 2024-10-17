from miRBench.encoder import get_encoder
from miRBench.predictor import get_predictor, list_predictors
from miRBench.dataset import get_dataset_path, list_datasets
import pandas as pd
import argparse
import os

def yield_blocks(file_path, block_size):
    with open(file_path, 'r') as file:
        # Read the first line as the header
        header = file.readline().strip().split('\t')
        block = []
        # Read from the second line onwards
        for line in file: 
            row = line.strip().split('\t')
            block.append(row)
            if len(block) == block_size:
                yield pd.DataFrame(block, columns=header)
                block = []
        if block:
            yield pd.DataFrame(block, columns=header)


def benchmark_all(block, dset, ratio):
    for tool in list_predictors():
        print(f"Running {tool} on {dset} dataset, ratio {ratio}")           
        encoder = get_encoder(tool)
        predictor = get_predictor(tool)
        input = encoder(block)
        output = predictor(input)
        block[tool] = output
    return block


def main():
    parser = argparse.ArgumentParser(description="Benchmark all available predictors on all available datasets")
    parser.add_argument("--out_dir", type=str, default=".", help="Output directory for predictions")
    
    args = parser.parse_args()
    
    split = "test"

    # loop over all available datasets
    for dset in list_datasets():
        for ratio in ["1", "10", "100"]:
            print(f"Downloading {dset} dataset, ratio {ratio}")
            input_file = get_dataset_path(dset, split=split, ratio=ratio)
            output_file = os.path.join(args.out_dir, f"{dset}_{ratio}_predictions.tsv")
            header_written = False
            for block in yield_blocks(input_file, 339066): 
                # 339066 is the number of examples in AGO2_eCLIP_Manakov2022_1_test_dataset.tsv, to be processed in one batch.
                # AGO2_eCLIP_Manakov2022_10_test_dataset.tsv and AGO2_eCLIP_Manakov2022_100_test_dataset.tsv will be processed in multiple batches.
                block_preds = benchmark_all(block, dset, ratio)
                if not header_written:
                    block_preds.to_csv(output_file, sep='\t', index=False, mode='w')
                    header_written = True
                else:
                    block_preds.to_csv(output_file, sep='\t', index=False, mode='a', header=False)

            print(f"Predictions for {dset} dataset, ratio {ratio}, written to {output_file}")

    print("Predictions for all datasets written to the output directory")

if __name__ == "__main__":
    main()