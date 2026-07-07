# workflow/rules/goflag_blast_locus_search.smk
# Summarize and filter BLAST results to identify sex-linked loci.

RESULTS_DIR = config["paths"]["results"]    # path for output files (result directory)
SAMPLES_PATH = config["paths"]["samples_tsv"]

checkpoint select_sex_linked_loci:
    """
    Analyze A_All reports from all samples to generate:
    1. A complete sex linked gene list (txt)
    2. Summary of shared sex linked loci by Order and species (tsv)
    """
    input:
        tsvs = expand(f"{RESULTS_DIR}/03_locus_search/{{sample_id}}/A_All_hits.tsv", sample_id=SAMPLES),
        samples_tsv = f"{SAMPLES_PATH}",
    output:
        out_list = f"{RESULTS_DIR}/04_filtered/all_sex_linked_genes.txt",
        out_order = f"{RESULTS_DIR}/04_filtered/sex_linked_summary_by_order.tsv",
        out_species = f"{RESULTS_DIR}/04_filtered/sex_linked_species_lists.tsv",
        out_autosomal = f"{RESULTS_DIR}/04_filtered/autosomal_genes.txt",
    log: "logs/4-1/sex_linked_loci_summary.log"
    shell:
        """
        python3 workflow/scripts/summarize_sex_linked.py \
            --tsvs {input.tsvs} \
            --samples_tsv {input.samples_tsv} \
            --out_list {output.out_list} \
            --out_order {output.out_order} \
            --out_species {output.out_species} \
            --out_autosomal {output.out_autosomal} \
            > {log} 2>&1
        """

# Run classify_duplications.py for each samples
rule classify_sample_duplications:
    """
    Reads the tagged BLAST hits for each sample and classifies multi-hit probes
    into specific duplication types (Local, WGD, Sex-linked, etc.).
    """
    input:
        tsv = f"{RESULTS_DIR}/03_locus_search/{{sample_id}}/A_All_hits.tsv"
    output:
        tsv = f"{RESULTS_DIR}/04_filtered/dup_tagged_og/{{sample_id}}_classified.tsv"
    params:
        max_dist = config["params"].get("max_dist", 100000)
    log:
        "logs/4-2/classify_duplications_{sample_id}.log"
    shell:
        """
        python3 workflow/scripts/classify_duplications.py \
            --input {input.tsv} \
            --output {output.tsv} \
            --max_dist {params.max_dist} \
            > {log} 2>&1
        """


# Create summary for all samples.
rule summarize_duplications_by_order:
    """
    Aggregates the classified duplication results from all samples, merges with
    metadata (samples.tsv) to summarize evolutionary events by taxonomic Order.
    """
    input:
        classified_tsvs = expand(f"{RESULTS_DIR}/04_filtered/dup_tagged_og/{{sample_id}}_classified.tsv", sample_id=SAMPLES),
        samples_tsv = f"{SAMPLES_PATH}",
    output:
        summary_by_order = f"{RESULTS_DIR}/04_filtered/duplication_summary_by_order.tsv",
        combined_all = f"{RESULTS_DIR}/04_filtered/all_samples_classified_duplications.tsv"
    log:
        "logs/4-3/summarize_duplications.log"
    shell:
        """
        python3 workflow/scripts/summarize_dup_type.py \
            --classified_tsvs {input.classified_tsvs} \
            --summary_by_order {output.summary_by_order} \
            --combined_all {output.combined_all} \
            > {log} 2>&1
        """
