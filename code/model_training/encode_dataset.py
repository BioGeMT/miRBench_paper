import pandas as pd
import numpy as np
import argparse
import time


def binding_encoding(df, alphabet={"AT": 1., "TA": 1., "GC": 1., "CG": 1.}, tensor_dim=(50, 20, 1),
                     ncRNA_col="noncodingRNA", gene_col="gene", label_col="label"):
    """
    Transform input sequence pairs to a binding matrix with corresponding labels.

    Parameters:
    - df: Pandas DataFrame with columns corresponding to ncRNA_col, gene_col, label_col
    - alphabet: dictionary with letter tuples as keys and 1s when they bind
    - tensor_dim: 2D binding matrix shape
    - ncRNA_col, gene_col, label_col: Column name for noncoding RNA sequences, gene sequences and label.

    Output:
    2D binding matrix, labels as np array
    """
    labels = df[label_col].to_numpy()

    # Initialize dot matrix with zeros
    ohe_matrix_2d = np.zeros((len(df), *tensor_dim), dtype="float32")

    df = df.reset_index(drop=True)

    # Compile matrix with Watson-Crick interactions
    for index, row in df.iterrows():
        for bind_index, bind_nt in enumerate(row[gene_col].upper()):
            if bind_index >= tensor_dim[0]:
                break
            for ncrna_index, ncrna_nt in enumerate(row[ncRNA_col].upper()):
                if ncrna_index >= tensor_dim[1]:
                    break
                base_pairs = bind_nt + ncrna_nt
                ohe_matrix_2d[index, bind_index, ncrna_index, 0] = alphabet.get(base_pairs, 0)
    return ohe_matrix_2d, labels


def encode_large_tsv_to_numpy(tsv_file_path, data_output_path, labels_output_path, chunk_size=10000,
                              ncRNA_col="noncodingRNA", gene_col="gene", label_col="label"):
    """
    Encode a large TSV file into a NumPy matrix using chunk processing.

    Parameters:
    - tsv_file_path: Path to the TSV file with dataset.
    - data_output_path: Path to the output data .npy file.
    - labels_output_path: Path to the output labels .npy file.
    - chunk_size: Number of rows to process at a time.
    - ncRNA_col, gene_col, label_col: Column name for noncoding RNA sequences, gene sequences and label.

    Output:
    The function just writes the encoded data and labels to the specified output files.
    """
    tensor_dim = (50, 20, 1)

    # Get total number of rows in the dataset
    num_rows = sum(len(df) for df in pd.read_csv(tsv_file_path, sep='\t', usecols=[0], chunksize=chunk_size))

    # Determine the shape of the output arrays
    labels_shape = (num_rows,)
    data_shape = (num_rows, *tensor_dim)

    try:
        # Create memory-mapped files
        ohe_matrix_2d = np.memmap(data_output_path, dtype='float32', mode='w+', shape=data_shape)
        labels = np.memmap(labels_output_path, dtype='float32', mode='w+', shape=labels_shape)

        row_offset = 0

        # Process each chunk
        for chunk in pd.read_csv(tsv_file_path, sep='\t', chunksize=chunk_size):
            encoded_data, encoded_labels = binding_encoding(chunk, ncRNA_col=ncRNA_col, gene_col=gene_col, label_col=label_col)

            # Write the chunk's data and labels to the memory-mapped files
            ohe_matrix_2d[row_offset:row_offset + len(chunk)] = encoded_data
            labels[row_offset:row_offset + len(chunk)] = encoded_labels
            row_offset += len(chunk)

        # Flush changes to disk
        ohe_matrix_2d.flush()
        labels.flush()
    except Exception as e:
        print(f"There was an unexpected error while encoding the dataset: {e}")
    finally:
        # Ensure the memory-mapped files are closed properly
        del ohe_matrix_2d
        del labels

def main():
    """
    Based on "miRBind: A deep learning method for miRNA binding classification." Genes 13.12 (2022): 2323. https://doi.org/10.3390/genes13122323.
    Original implementation: https://github.com/ML-Bioinfo-CEITEC/miRBind

    Encodes miRNA and gene sequences into 2D-binding matrix.
    2D-binding matrix has shape (gene_max_len=50, miRNA_max_len=20, 1) and contains 1 for Watson-Crick interactions and 0 otherwise.
    """

    parser = argparse.ArgumentParser(
        description="Encode dataset to miRNA x target binding matrix. Outputs numpy file with matrices and and numpy file with corresponding labels. Expected columns of the dataset are 'noncodingRNA', 'gene' and 'label'")
    parser.add_argument('-i', '--i_file', type=str, required=True, help="Input dataset file name")
    parser.add_argument('-o', '--o_prefix', type=str, required=True, help="Output file name prefix")
    parser.add_argument('--ncRNA_column', type=str, default='noncodingRNA', help="Name of the column with noncoding RNA sequences")
    parser.add_argument('--gene_column', type=str, default='gene', help="Name of the column with gene sequences")
    parser.add_argument('--label_column', type=str, default='label', help="Name of the column with labels")

    args = parser.parse_args()

    start = time.time()
    encode_large_tsv_to_numpy(args.i_file, args.o_prefix + '_dataset.npy', args.o_prefix + '_labels.npy',
                              ncRNA_col=args.ncRNA_column, gene_col=args.gene_column, label_col=args.label_column)
    end = time.time()

    print("Elapsed time: ", end - start, " s.")


if __name__ == "__main__":
    main()