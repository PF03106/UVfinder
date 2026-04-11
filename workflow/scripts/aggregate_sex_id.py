#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Aggregate sex assignment results and merge with metadata.")
    
    parser.add_argument("--inputs", required=True, nargs='+', help="List of input TSV files from 02_assign_sex.smk")
    parser.add_argument("--metadata", required=True, help="Path to the samples.tsv metadata file")
    parser.add_argument("--output", required=True, help="Path for the final aggregated output TSV file")
    
    args = parser.parse_args()

    # 1. combine all .tsv files into a single DataFrame
    dfs = [pd.read_csv(f, sep='\t') for f in args.inputs]
    combined_df = pd.concat(dfs, ignore_index=True)

    meta_df = pd.read_csv(args.metadata, sep='\t')
    final_df = pd.merge(combined_df, meta_df[['sample_id', 'genus', 'species']], on='sample_id', how='left')
    cols = ['sample_id', 'genus', 'species'] + [col for col in final_df.columns if col not in ['sample_id', 'genus', 'species']]
    final_df = final_df[cols]
    
    # 2. save the final aggregated DataFrame to a TSV file
    final_df.to_csv(args.output, sep='\t', index=False)
    
    print(f"Successfully aggregated {len(args.inputs)} files into {args.output}")

if __name__ == "__main__":
    main()