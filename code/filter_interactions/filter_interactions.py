import pandas as pd
import argparse

def filter_interactions(df):
    # Canonical seed: Seed6mer is 1    
    df_canonical = df[df['Seed6mer'] == 1].copy()
    df_canonical = df_canonical.drop(columns=['Seed6mer', 'Seed6merBulgeOrMismatch'])

    # Note that in miRBench package, Seed6merBulgeOrMismatch is inclusive of Seed6mer
    # Non-canonical seed: Seed6merBulgeOrMismatch is 1 AND Seed6mer is 0
    df_noncanonical = df.loc[(df["Seed6merBulgeOrMismatch"] == 1) & (df["Seed6mer"] == 0)].copy()
    df_noncanonical = df_noncanonical.drop(columns=['Seed6mer', 'Seed6merBulgeOrMismatch'])

    # No seed: Seed6merBulgeOrMismatch is 0
    df_noseed = df[df['Seed6merBulgeOrMismatch'] == 0].copy()
    df_noseed = df_noseed.drop(columns=['Seed6mer', 'Seed6merBulgeOrMismatch'])

    return df_canonical, df_noncanonical, df_noseed

def write_interactions(df, ofile):
    df.to_csv(ofile, sep='\t', index=False)
        
def main():
    parser = argparse.ArgumentParser(description="Filter canonical/non-canonical/no-seed interactions, for all Manakov datasets")
    parser.add_argument("--ifile", type=str, help="Input file with seed types")
    parser.add_argument("--canonical_ofile", type=str, help="Output file for canonical seed types")
    parser.add_argument("--noncanonical_ofile", type=str, help="Output file for noncanonical seed types")
    parser.add_argument("--nonseed_ofile", type=str, help="Output file for nonseed types")
    args = parser.parse_args()

    # Read file with seed types
    df_seed_types = pd.read_csv(args.ifile, sep='\t')
    
    # Filter canonical/non-canonical/non-seed interactions
    df_canonical, df_noncanonical, df_noseed = filter_interactions(df_seed_types)

    # Write interactions to file
    write_interactions(df_canonical, args.canonical_ofile)
    write_interactions(df_noncanonical, args.noncanonical_ofile)
    write_interactions(df_noseed, args.nonseed_ofile)

if __name__ == "__main__":
    main()