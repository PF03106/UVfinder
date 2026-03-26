# workflow/rules/goflag_blast_locus_search.smk
# Summarize and filter BLAST results to identify sex-linked loci.

RESULTS_DIR = config["paths"]["results"]    # path for output files (result directory)

checkpoint select_sex_linked_loci:
    """
    Analyze A_All reports from all samples to generate:
    1. A complete sex linked gene list (txt)
    2. Summary of shared sex linked loci by Order and species (tsv)
    """
    input:
        tsvs = expand(f"{RESULTS_DIR}/03_locus_search/{{sample_id}}/A_All_hits.tsv", sample_id=SAMPLES),
        samples_tsv = "config/samples.tsv"
    output:
        out_list = f"{RESULTS_DIR}/04_filtered/all_sex_linked_genes.txt",
        out_order = f"{RESULTS_DIR}/04_filtered/sex_linked_summary_by_order.tsv",
        out_species = f"{RESULTS_DIR}/04_filtered/sex_linked_species_lists.tsv"
    log: "logs/4/sex_linked_loci_summary.log"
    shell:
        """
        python3 workflow/scripts/summarize_sex_linked.py \
            --tsvs {input.tsvs} \
            --samples_tsv {input.samples_tsv} \
            --out_list {output.out_list} \
            --out_order {output.out_order} \
            --out_species {output.out_species} \
            > {log} 2>&1
        """