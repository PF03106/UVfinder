# Phase 4: Locus Filtering (Expanded Census using All Hits)

rule select_interesting_loci:
    """
    Step 4.1: Perform a full census across ALL hits.
    A locus is selected if ANY hit in the A_all folder is tagged with _U or _V.
    """
    input:
        # Adjust A or B
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

rule collect_and_qc_loci:
    """
    Step 4.2: Collect all sequences from A_all for the selected loci and perform QC.
    This captures all duplicated copies for detailed evolutionary analysis.
    """
    input:
        loci_list = rules.select_interesting_loci.output.list,
        all_hits = expand("results/03_locus_search/{s}/A_all", s=SAMPLES)
    output:
        qc_passed = directory("results/04_filtered/qc_passed_loci")
    log: "logs/phase4_qc.log"
    shell:
        """
        python3 workflow/scripts/qc_sequences.py \
            --loci_list {input.loci_list} \
            --input_dirs {input.all_hits} \
            --output_dir {output.qc_passed}
        """