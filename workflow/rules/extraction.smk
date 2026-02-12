# Phase 5: Extract Sequences (Best & All)

# 1. Extract Best Hits
rule extract_best_sequences:
    input:
        tsv = "results/03_locus_search/{s}/B_Best_hits.tsv",
        genome = "results/00_renamed/{s}_renamed.fasta",
        loci_list = "results/04_filtered/all_sex_linked_genes.txt"
    output:
        out_dir = directory("results/05_extracted/{s}/Best")
    params:
        sample = "{s}",
        # Get values from config (use default if not set)
        flank = config.get("flanking_bp", 20),
        # Best hits only have 1 record, but for code compatibility we use max_hits parameter
        max_hits = config.get("max_blast_hits", 10)
    log: "logs/5_best/extract_best_{s}.log"
    shell:
        """
        python3 workflow/scripts/extract_from_tsv.py \
            --tsv {input.tsv} \
            --genome {input.genome} \
            --loci_list {input.loci_list} \
            --out_dir {output.out_dir} \
            --sample {params.sample} \
            --flank {params.flank} \
            --max_hits {params.max_hits} \
            > {log} 2>&1
        """

# 2. Extract All Hits
rule extract_all_sequences:
    input:
        tsv = "results/03_locus_search/{s}/A_All_hits.tsv",
        genome = "results/00_renamed/{s}_renamed.fasta",
        loci_list = "results/04_filtered/all_sex_linked_genes.txt"
    output:
        out_dir = directory("results/05_extracted/{s}/All")
    params:
        sample = "{s}",
        flank = config.get("flanking_bp", 20),
        max_hits = config.get("max_blast_hits", 10)
    log: "logs/5_all/extract_all_{s}.log"
    shell:
        """
        python3 workflow/scripts/extract_from_tsv.py \
            --tsv {input.tsv} \
            --genome {input.genome} \
            --loci_list {input.loci_list} \
            --out_dir {output.out_dir} \
            --sample {params.sample} \
            --flank {params.flank} \
            --max_hits {params.max_hits} \
            > {log} 2>&1
        """