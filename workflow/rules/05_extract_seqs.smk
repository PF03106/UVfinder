# workflow/rules/05_extract_seqs.smk
# Extract sequences for goflag probes from BLAST results with flanking regions set by users.

RESULTS_DIR = config["paths"]["results"]    # path for output files (result directory)

# 1. Extract Best Hits
rule extract_best_sequences:
    input:
        tsv = f"{RESULTS_DIR}/03_locus_search/{{sample_id}}/B_Best_hits.tsv",
        genome = f"{RESULTS_DIR}/00_renamed/{{sample_id}}_renamed.fasta",
        loci_list = f"{RESULTS_DIR}/04_filtered/all_sex_linked_genes.txt"
    output:
        out_dir = directory(f"{RESULTS_DIR}/05_extracted/{{sample_id}}/Best")
    params:
        sample = "{sample_id}",
        flank = config["params"]["flanking_bp"],
        # Best hits only have 1 record, but for code compatibility we use max_hits parameter
        max_hits = config["params"]["max_blast_hits"]
    log: "logs/5_best/extract_best_{sample_id}.log"
    shell:
        """
        python3 workflow/scripts/extract_sex_linked_with_flanks.py \
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
        tsv = f"{RESULTS_DIR}/03_locus_search/{{sample_id}}/A_All_hits.tsv",
        genome = f"{RESULTS_DIR}/00_renamed/{{sample_id}}_renamed.fasta",
        loci_list = f"{RESULTS_DIR}/04_filtered/all_sex_linked_genes.txt"
    output:
        out_dir = directory(f"{RESULTS_DIR}/05_extracted/{{sample_id}}/All")
    params:
        sample = "{sample_id}",
        flank = config["params"]["flanking_bp"],
        max_hits = config["params"]["max_blast_hits"]
    log: "logs/5_all/extract_all_{sample_id}.log"
    shell:
        """
        python3 workflow/scripts/extract_sex_linked_with_flanks.py \
            --tsv {input.tsv} \
            --genome {input.genome} \
            --loci_list {input.loci_list} \
            --out_dir {output.out_dir} \
            --sample {params.sample} \
            --flank {params.flank} \
            --max_hits {params.max_hits} \
            > {log} 2>&1
        """

# 3. Extract Autosomal Hits (Non-target, genes that are not in all_sex_linked_genes.txt)
rule extract_not_sex_linked_sequences:
    input:
        tsv = f"{RESULTS_DIR}/03_locus_search/{{sample_id}}/B_Best_hits.tsv",
        genome = f"{RESULTS_DIR}/00_renamed/{{sample_id}}_renamed.fasta",
        loci_list = f"{RESULTS_DIR}/04_filtered/all_sex_linked_genes.txt"
    output:
        out_dir = directory(f"{RESULTS_DIR}/05_extracted/{{sample_id}}/Not_sex_linked")
    params:
        sample = "{sample_id}",
        flank = config["params"]["flanking_bp"],
        max_hits = config["params"]["max_blast_hits"]
    log: "logs/5_best/extract_best_nontarget_{sample_id}.log"
    shell:
        """
        python3 workflow/scripts/extract_sex_linked_with_flanks.py \
            --tsv {input.tsv} \
            --genome {input.genome} \
            --loci_list {input.loci_list} \
            --out_dir {output.out_dir} \
            --sample {params.sample} \
            --flank {params.flank} \
            --max_hits {params.max_hits} \
            --exclude \
            > {log} 2>&1
        """