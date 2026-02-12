# workflow/rules/goflag_blast_locus_search.smk

# 1. Call Python script to generate order-specific probes
rule prepare_order_probes:
    input:
        probe_dir = "resources/query_sets", 
        samples_tsv = "config/samples.tsv"
    output:
        # Sample-specific probe file
        # Using temp() will automatically delete after BLAST completes to save space (optional)
        #specific_probes = temp("resources/probes_by_sample/{s}_specific_probes.fasta")
        specific_probes = "resources/probes_by_sample/{s}_specific_probes.fasta"
    log: "logs/3-1/prepare_probes_{s}.log"
    shell:
        """
        python3 workflow/scripts/blast_by_orders.py \
            --sample {wildcards.s} \
            --samples_tsv {input.samples_tsv} \
            --probe_dir {input.probe_dir} \
            --output {output.specific_probes} \
            > {log} 2>&1
        """

# 2. Run BLAST with custom probes
rule blast_search:
    input:
        genome = "results/00_renamed/{s}_renamed.fasta",
        # Input file created above
        query = "resources/probes_by_sample/{s}_specific_probes.fasta"
    output:
        blast_out = "results/03_locus_search/{s}/{s}_blast_results.txt"
    threads: 8
    params:
        evalue = config["params"]["query_blast_evalue"],
        outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    log: "logs/3-2/blast_{s}.log"
    shell:
        """
        tblastn \
            -query {input.query} \
            -subject {input.genome} \
            -evalue {params.evalue} \
            -outfmt "{params.outfmt}" \
            -num_threads {threads} \
            -out {output.blast_out} \
            2> {log}
        """

# 3. Partition BLAST results into All hits and Best hits. 
rule partition_blast_results:
    """
    Step 3.3: Process raw BLAST results to create organized tables (TSV).
    - Maps BLAST hits to Sex Chromosomes (U, V, Autosome, Unknown).
    - Ranks hits by E-value (Rank 1 = Best).
    - Splits results into 'Best' (Rank 1 only) and 'All' (All ranks).
    * Note: This step does NOT extract sequences yet.
    """
    input:
        blast = "results/03_locus_search/{s}/{s}_blast_results.txt",
        sex_map = "results/02_sex_id/{s}_sex_assignment.tsv"
    output:
        best_tsv = "results/03_locus_search/{s}/B_Best_hits.tsv",
        all_tsv = "results/03_locus_search/{s}/A_All_hits.tsv"
    log: "logs/3-3/All_or_Best_partition_{s}.log"
    shell:
        """
        python3 workflow/scripts/extract_and_tag.py \
            --blast {input.blast} \
            --sex_map {input.sex_map} \
            --out_best {output.best_tsv} \
            --out_all {output.all_tsv} \
            > {log} 2>&1
        """
