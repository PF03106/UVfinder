# workflow/rules/06_align_trim.smk
# Phase 6: Multiple Sequence Alignment (MSA) using MAFFT add & Trimming using TrimAl.

RESULTS_DIR = config["paths"]["results"]    # path for output files (result directory)

# Rule 1: Gather Sequences for MSA (Merge sequences from all samples for each gene)
rule gather_sequences:
    input:
        extracted_dirs = lambda wildcards: expand(f"{RESULTS_DIR}/05_extracted/{{sample_id}}/{wildcards.type}", sample_id=SAMPLES)
    output:
        merged = f"{RESULTS_DIR}/06_alignment/{{type}}/merged/{{gene}}.fasta"
    params:
        base_dir = f"{RESULTS_DIR}/05_extracted",
        samples = lambda w: " ".join(SAMPLES),
        min_taxa = config["params"]["mafft"]["min_taxa"]
    log: "logs/6-1/{type}/gather_{gene}.log"
    shell:
        """
        python3 workflow/scripts/gather_seqs_for_msa.py \
            --base_dir {params.base_dir} \
            --out_file {output.merged} \
            --samples {params.samples} \
            --type_dir {wildcards.type} \
            --gene_id {wildcards.gene} \
            --min_taxa {params.min_taxa} \
            > {log} 2>&1
        """
# Rule 2: MAFFT Alignment (Alignment with Reference included)
rule mafft_alignment:
    input:
        merged_fasta = f"{RESULTS_DIR}/06_alignment/{{type}}/merged/{{gene}}.fasta"
    params:
        ref_path = lambda w: os.path.join(config["paths"]["ref_alignment_dir"], f"{w.gene}.fasta"),
        mafft_opts = config["params"]["mafft"]["alignment_option"]
    output:
        aln_ref = temp(f"{RESULTS_DIR}/06_alignment/{{type}}/aligned/{{gene}}_with_ref.aln")
    threads: 4
    log: "logs/6-2/mafft_{type}_{gene}.log"
    shell:
        """
        if [ -s "{input.merged_fasta}" ]; then
            if [ -f "{params.ref_path}" ]; then
                echo "Reference found: {params.ref_path}" > {log}
                # Add to Reference (--add)
                mafft {params.mafft_opts} --add {input.merged_fasta} \
                      --reorder {params.ref_path} > {output.aln_ref} 2>> {log}
            else
                echo "No reference found. Running de novo." >> {log}
                mafft {params.mafft_opts} \
                      {input.merged_fasta} > {output.aln_ref} 2>> {log}
            fi
        else
            touch {output.aln_ref}
        fi
        """

# Rule 3: Filter Alignment
rule filter_alignment:
    input:
        # MAFFT result (Ref included)
        aln_ref = f"{RESULTS_DIR}/06_alignment/{{type}}/aligned/{{gene}}_with_ref.aln"
    output:
        # Final alignment file (Ref removed)
        aln_final = temp(f"{RESULTS_DIR}/06_alignment/{{type}}/temp_aligned/{{gene}}.aln")
    params:
        # Pass my sample list
        samples = lambda w: " ".join(SAMPLES)
    log: "logs/6-3/filter_{type}_{gene}.log"
    shell:
        """
        if [ -s "{input.aln_ref}" ]; then
            python3 workflow/scripts/filter_alignments.py \
                {input.aln_ref} \
                {output.aln_final} \
                {params.samples} \
                > {log} 2>&1
        else
            touch {output.aln_final}
        fi
        """

# Rule 4: TrimAl (Trimming with file containing only my samples)
rule trimal_trimming:
    input:
        # Input file after filtering
        aln = f"{RESULTS_DIR}/06_alignment/{{type}}/temp_aligned/{{gene}}.aln"
    output:
        trimmed = f"{RESULTS_DIR}/06_alignment/{{type}}/trimmed/{{gene}}.trimmed.aln"
    params:
        trimal_opts = config["params"]["trimal"]
    log: "logs/6-4/trimal_{type}_{gene}.log"
    shell:
        """
        if [ -s "{input.aln}" ]; then
            trimal -in {input.aln} -out {output.trimmed} {params.trimal_opts} > {log} 2>&1 || true
            if [ ! -f "{output.trimmed}" ]; then
                echo "trimAl produced no output for {wildcards.gene}. Creating empty file." >> {log}
                touch {output.trimmed}
            fi
        else
            echo "Input alignment is empty. Creating empty output." > {log}
            touch {output.trimmed}
        fi
        """