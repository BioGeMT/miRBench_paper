import argparse
import sys
import warnings

import pandas as pd


INPUT_COLUMNS = [
    "gene",
    "noncodingRNA",
    "noncodingRNA_name",
    "noncodingRNA_fam",
    "feature",
    "test",
    "label",
    "chr",
    "start",
    "end",
    "strand",
    "Nunique",
    "dominant_region_selected",
    "regions_present_selected",
    "selected_read_start_in_tx_1based",
    "selected_read_end_in_tx_1based",
]

OUTPUT_COLUMNS = [
    "gene",
    "noncodingRNA",
    "noncodingRNA_name",
    "noncodingRNA_fam",
    "feature",
    "test",
    "label",
    "chr",
    "start",
    "end",
    "strand",
    "Nunique",
    "dominant_region",
    "regions_present",
    "read_start_in_sel_tx_1based",
    "read_end_in_sel_tx_1based",
]

RENAME_MAP = {
    "dominant_region_selected": "dominant_region",
    "regions_present_selected": "regions_present",
    "selected_read_start_in_tx_1based": "read_start_in_sel_tx_1based",
    "selected_read_end_in_tx_1based": "read_end_in_sel_tx_1based",
}

INT_COLUMNS = [
    "read_start_in_sel_tx_1based",
    "read_end_in_sel_tx_1based",
]


def read_input(input_file):
    try:
        if input_file:
            return pd.read_csv(input_file, sep="\t")
        return pd.read_csv(sys.stdin, sep="\t")
    except Exception as exc:
        warnings.warn(f"Error reading input data: {exc}", category=UserWarning)
        sys.exit(1)


def filter_columns(data):
    missing_columns = [column for column in INPUT_COLUMNS if column not in data.columns]
    if missing_columns:
        warnings.warn(
            "Missing required columns: " + ", ".join(missing_columns),
            category=UserWarning,
        )
        sys.exit(1)

    filtered_data = data.loc[:, INPUT_COLUMNS].rename(columns=RENAME_MAP)

    for column in INT_COLUMNS:
        filtered_data[column] = pd.to_numeric(filtered_data[column], errors="coerce").astype(
            "Int64"
        )

    return filtered_data.loc[:, OUTPUT_COLUMNS]


def write_output(data, output_file):
    try:
        if output_file:
            data.to_csv(output_file, sep="\t", index=False)
        else:
            data.to_csv(sys.stdout, sep="\t", index=False)
    except Exception as exc:
        warnings.warn(f"Error writing output data: {exc}", category=UserWarning)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Filter genomic region annotator site summaries to the training columns."
    )
    parser.add_argument("--ifile", help="Input file (default: STDIN)")
    parser.add_argument("--ofile", help="Output file (default: STDOUT)")
    args = parser.parse_args()

    data = read_input(args.ifile)
    filtered_data = filter_columns(data)
    write_output(filtered_data, args.ofile)


if __name__ == "__main__":
    main()
