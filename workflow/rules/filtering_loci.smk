# Phase 4: Summarize & Filter Sex-linked Loci

checkpoint select_sex_linked_loci:
    """
    Analyze A_All reports from all samples to generate:
    1. A complete sex linked gene list (txt)
    2. Summary of shared sex linked loci by Order (tsv)
    """
    input:
        tsvs = expand("results/03_locus_search/{s}/A_All_hits.tsv", s=SAMPLES),
        samples_tsv = "config/samples.tsv"
    output:
        out_list = "results/04_filtered/all_sex_linked_genes.txt",
        out_order = "results/04_filtered/sex_linked_summary_by_order.tsv",
        out_species = "results/04_filtered/sex_linked_species_lists.tsv"
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