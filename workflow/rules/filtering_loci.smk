# Phase 4: Locus Filtering & Collection (No QC)

checkpoint select_interesting_loci:
    """
    Step 4.1: Checkpoint rule to determine dynamic outputs.
    """
    input:
        all_hits = expand("results/03_locus_search/{s}/A_all", s=SAMPLES)
    output:
        list = "results/04_filtered/interesting_loci.txt"
    log: "logs/phase4_select_loci.log"
    shell:
        """
        python3 workflow/scripts/select_loci.py \
            --input_dirs {input.all_hits} \
            --output {output.list}
        """
        
rule collect_interesting_loci:
    """
    Step 4.2: Collect all sequences for selected loci (WITHOUT QC).
    """
    input:
        loci_list = rules.select_interesting_loci.output.list,
        all_hits = expand("results/03_locus_search/{s}/A_all", s=SAMPLES)
    output:
        collection = directory("results/04_filtered/collected_loci")
    log: "logs/phase4_collection.log"
    shell:
        """
        python3 workflow/scripts/collect_sequences.py \
            --loci_list {input.loci_list} \
            --input_dirs {input.all_hits} \
            --output_dir {output.collection}
        """