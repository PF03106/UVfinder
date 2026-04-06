#!/usr/bin/env python3
# Phase 4: Summarize & Filter Sex-linked Loci
import pandas as pd
import argparse
import os

def summarize_sex_linked(tsv_files, samples_tsv, out_list, out_not_sex_linked, out_order, out_species):
    print(f"--- Summarizing Sex-linked Loci from {len(tsv_files)} samples ---")
    
    # 1. Load sample metadata (Order)
    try:
        sample_meta = pd.read_csv(samples_tsv, sep='\t')
        sample_info = sample_meta.set_index('sample_id')[['order', 'genus', 'species']].to_dict('index')
    except Exception as e:
        print(f"❌ Error reading samples.tsv: {e}")
        return

    collected_data_dict = {}
    all_genes = set()

    # 2. Iterate through all sample TSVs
    for tsv in tsv_files:
        if not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
            continue
            
        sample_id = os.path.basename(os.path.dirname(tsv))
        info = sample_info.get(sample_id, {'order': 'Unknown', 'genus': 'Unknown', 'species': 'Unknown'})
        order = info['order']
        genus = info['genus']
        species = info['species']

        try:
            df = pd.read_csv(tsv, sep='\t')
            if 'sex_tag' not in df.columns: continue

            if 'qseqid' in df.columns:
                current_genes = df['qseqid'].astype(str).apply(lambda x: x.split('|')[0])
                all_genes.update(current_genes)

            targets = df[df['sex_tag'].str.upper().isin(['U', 'V'])]
            
            if not targets.empty:
                for _, row in targets.iterrows():
                    locus = str(row['qseqid']).split('|')[0]
                    sex_tag = row['sex_tag']
                    
                    key = (locus, sample_id, order, sex_tag)
                    collected_data_dict[key] = {
                        "Gene": locus,
                        "Sample": sample_id,
                        "Order": order,
                        "Genus": genus,
                        "Species": species,
                        "Sex_Tag": sex_tag
                    }
                    
        except Exception as e:
            print(f"⚠️ Error processing {tsv}: {e}")

    # If no data, create empty files
    if not collected_data_dict:
        open(out_list, 'w').close()
        open(out_not_sex_linked, 'w').close()
        pd.DataFrame(columns=["Order", "Gene", "Count", "Samples"]).to_csv(out_order, sep='\t', index=False)
        pd.DataFrame(columns=["Sample", "Order","Genus", "Species", "Gene_Count", "Sex_Linked_Genes"]).to_csv(out_species, sep='\t', index=False)
        return

    full_df = pd.DataFrame(list(collected_data_dict.values()))

    # --- [File 1] Gene list (all_sex_linked_genes.txt) ---
    unique_genes = sorted(full_df["Gene"].unique())
    with open(out_list, 'w') as f:
        for gene in unique_genes:
            f.write(f"{gene}\n")

    # --- [File 4] Not Sex-linked Gene list (not_sex_linked_genes.txt) ---
    not_sex_linked_genes = sorted(list(all_genes - set(unique_genes)))
    with open(out_not_sex_linked, 'w') as f:
        for gene in not_sex_linked_genes:
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
    pd.DataFrame(summary_list).sort_values(by=["Order", "Count"], ascending=[True, False]).to_csv(out_order, sep='\t', index=False)

    # --- [File 3] Single TSV summary by Sample/Species ---
    sample_summary_list = []
    for sample_name, group in full_df.groupby("Sample"):
        genes = sorted(group["Gene"].unique())
        
        sample_summary_list.append({
            "Sample": sample_name,
            "Gene_Count": len(genes),
            "Sex_Linked_Genes": ",".join(genes)
        })
        
    sample_summary_df = pd.DataFrame(sample_summary_list)
    all_samples_meta = sample_meta[['sample_id', 'order', 'genus', 'species']].rename(columns={'sample_id': 'Sample', 'order': 'Order', 'genus': 'Genus', 'species': 'Species'})
    final_sample_summary = pd.merge(all_samples_meta, sample_summary_df, on='Sample', how='left')
    final_sample_summary['Gene_Count'] = final_sample_summary['Gene_Count'].fillna(0).astype(int)
    final_sample_summary['Sex_Linked_Genes'] = final_sample_summary['Sex_Linked_Genes'].fillna("")
    columns_order = ["Sample", "Order", "Genus", "Species", "Gene_Count", "Sex_Linked_Genes"]
    final_sample_summary[columns_order].to_csv(out_species, sep='\t', index=False)

    print(f"✅ Summary Generated:")
    print(f"   - Total Unique Genes across all hits: {len(all_genes)}")
    print(f"   - Sex-linked Genes: {len(unique_genes)}")
    print(f"   - Not Sex-linked Genes: {len(not_sex_linked_genes)}")
    print(f"   - Output 1 (Sex-linked List): {out_list}")
    print(f"   - Output 2 (By Order): {out_order}")
    print(f"   - Output 3 (By Species): {out_species}")
    print(f"   - Output 4 (Not Sex-linked List): {out_not_sex_linked}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsvs", nargs='+', required=True)
    parser.add_argument("--samples_tsv", required=True)
    parser.add_argument("--out_list", required=True)
    parser.add_argument("--out_not_sex_linked", required=True)
    parser.add_argument("--out_order", required=True)
    parser.add_argument("--out_species", required=True) 
    
    args = parser.parse_args()
    
    summarize_sex_linked(args.tsvs, args.samples_tsv, args.out_list, args.out_not_sex_linked, args.out_order, args.out_species)