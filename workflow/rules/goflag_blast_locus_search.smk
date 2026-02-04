import os
from Bio import SeqIO

# --- Phase 3: Whole-Genome Locus Search (GoFlag Probes) ---
# This phase searches for ~400 orthologous loci across the genomes,
# tags them based on sex-linked status, and extracts sequences.

rule prepare_sample_queries:
    """
    Step 3.1: Extract probe sequences that match the Order of the current sample.
    This ensures we use evolutionarily close sequences for BLAST, increasing sensitivity.
    """
    input:
        probe_dir = "resources/query_sets" # Folder containing the 400 refined FASTA files
    output:
        temp_fasta = temp("results/03_locus_search/{sample_id}/matched_queries.fasta")
    run:
        # Retrieve the Order of the current sample from the samples_df dataframe.
        # We use the dataframe lookup because the filename might be short (e.g., S0046).
        try:
            # samples_df is inherited from the main Snakefile
            sample_order = samples_df.loc[wildcards.sample_id, "order"]
        except KeyError:
            print(f"❌ Error: Sample ID '{wildcards.sample_id}' not found in samples.tsv index.")
            raise

        # Normalize the order name (lowercase and remove spaces) to match refined probe headers.
        order_norm = str(sample_order).lower().replace(" ", "").strip()
        
        print(f"--- 🧬 Extracting GoFlag probes for {wildcards.sample_id} (Target Order: {order_norm}) ---")
        
        count = 0
        with open(output.temp_fasta, "w") as out_f:
            # Iterate through the 400 refined probe FASTA files
            for locus_file in os.listdir(input.probe_dir):
                if not locus_file.endswith(".fasta"):
                    continue
                
                locus_path = os.path.join(input.probe_dir, locus_file)
                for record in SeqIO.parse(locus_path, "fasta"):
                    # Check if the probe header (format: >order_species) matches our sample's order
                    if record.id.lower().startswith(f"{order_norm}_"):
                        SeqIO.write(record, out_f, "fasta")
                        count += 1
        
        if count == 0:
            print(f"⚠️ Warning: No matching probes found for order '{order_norm}' in {wildcards.sample_id}.")

rule blast_queries:
    """
    Step 3.2: Perform tBLASTn using the order-specific queries against the sample genome.
    """
    input:
        query = rules.prepare_sample_queries.output.temp_fasta,
        db_index = "results/01_blastdb/{sample_id}_renamed.nhr"
    output:
        "results/03_locus_search/{sample_id}_goflag.tblastn"
    params:
        db_prefix = "results/01_blastdb/{sample_id}_renamed",
        evalue = config["params"]["query_blast_evalue"] # Threshold defined in config.yaml
    shell:
        """
        tblastn -query {input.query} \
                -db {params.db_prefix} \
                -evalue {params.evalue} \
                -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
                -out {output}
        """

rule filter_and_tag_hits:
    """
    Step 3.3: Assign sex-linked tags (_U-linked, _V-linked, or _A) based on contig mapping.
    Outputs are separated into 'All hits' (A) and 'Best hits' (B) folders.
    """
    input:
        blast_res = rules.blast_queries.output,
        sex_map = "results/02_sex_id/{sample_id}_sex_assignment.tsv", # Identification results from Phase 2
        genome = "results/00_renamed/{sample_id}_renamed.fasta"
    output:
        all_hits = directory("results/03_locus_search/{sample_id}/A_all"),
        best_hits = directory("results/03_locus_search/{sample_id}/B_Best")
    shell:
        """
        python3 workflow/scripts/extract_and_tag.py \
            --blast {input.blast_res} \
            --sex_map {input.sex_map} \
            --genome {input.genome} \
            --out_best {output.best_hits} \
            --out_all {output.all_hits}
        """