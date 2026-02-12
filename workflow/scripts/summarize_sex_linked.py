import pandas as pd
import argparse
import os

def summarize_sex_linked(tsv_files, samples_tsv, out_list, out_order):
    print(f"--- Summarizing Sex-linked Loci from {len(tsv_files)} samples ---")
    
    # 1. Load sample metadata (Order)
    try:
        sample_meta = pd.read_csv(samples_tsv, sep='\t')
        sample_to_order = sample_meta.set_index('sample_id')['order'].to_dict()
    except Exception as e:
        print(f"❌ Error reading samples.tsv: {e}")
        return

    # Dictionary for efficient duplicate removal (key: (Gene, Sample, Order, Sex_Tag))
    collected_data_dict = {}

    # 2. Iterate through all sample TSVs
    for tsv in tsv_files:
        if not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
            continue
            
        sample_id = os.path.basename(os.path.dirname(tsv))
        order = sample_to_order.get(sample_id, "Unknown")
        
        try:
            df = pd.read_csv(tsv, sep='\t')
            if 'sex_tag' not in df.columns: continue

            # Filter for U/V sex tags
            targets = df[df['sex_tag'].str.upper().isin(['U', 'V'])]
            
            if not targets.empty:
                for _, row in targets.iterrows():
                    # Extract gene ID cleanly (G4989|... -> G4989)
                    locus = str(row['qseqid']).split('|')[0]
                    sex_tag = row['sex_tag']
                    
                    # Use tuple as unique key to prevent duplicates at insertion time
                    key = (locus, sample_id, order, sex_tag)
                    collected_data_dict[key] = {
                        "Gene": locus,
                        "Sample": sample_id,
                        "Order": order,
                        "Sex_Tag": sex_tag
                    }
                    
        except Exception as e:
            print(f"⚠️ Error processing {tsv}: {e}")

    # If no data, create empty files
    if not collected_data_dict:
        open(out_list, 'w').close()
        pd.DataFrame(columns=["Order", "Gene", "Count", "Samples"]).to_csv(out_order, sep='\t', index=False)
        return

    # Convert dictionary values to DataFrame (already de-duplicated)
    full_df = pd.DataFrame(list(collected_data_dict.values()))

    # --- [File 1] Gene list (all_sex_linked_genes.txt) ---
    unique_genes = sorted(full_df["Gene"].unique())
    with open(out_list, 'w') as f:
        for gene in unique_genes:
            f.write(f"{gene}\n")
            
    # --- [File 2] Statistical summary by Order (sex_linked_summary_by_order.tsv) ---
    summary_list = []
    
    for (order_name, gene_name), group in full_df.groupby(["Order", "Gene"]):
        count = group["Sample"].nunique()
        samples_str = ",".join(sorted(group["Sample"].unique()))
        
        summary_list.append({
            "Order": order_name,
            "Gene": gene_name,
            "Count": count,
            "Samples": samples_str
        })
    
    summary_df = pd.DataFrame(summary_list)
    summary_df = summary_df.sort_values(by=["Order", "Count"], ascending=[True, False])
    summary_df.to_csv(out_order, sep='\t', index=False)

    print(f"✅ Summary Generated:")
    print(f"   - Unique Genes: {len(unique_genes)}")
    print(f"   - Output: {out_order}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsvs", nargs='+', required=True)
    parser.add_argument("--samples_tsv", required=True)
    parser.add_argument("--out_list", required=True)
    parser.add_argument("--out_order", required=True)
    
    args = parser.parse_args()
    
    summarize_sex_linked(args.tsvs, args.samples_tsv, args.out_list, args.out_order)