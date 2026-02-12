# Phase 6: Multiple Sequence Alignment (MSA) & Trimming

# Rule 6-1: Gather Sequences
rule gather_sequences:
    input:
        extracted_dirs = expand("results/05_extracted/{s}/{type}", s=SAMPLES, type=["Best", "All"]),
        gene_list = "results/04_filtered/all_sex_linked_genes.txt"
    output:
        merged = "results/06_alignment/{type}/merged/{gene}.fasta"
    params:
        base_dir = "results/05_extracted",
        samples = lambda w: " ".join(SAMPLES),
        min_taxa = config["params"]["mafft"]["min_taxa"]
    log: "logs/6-1/gather_{type}_{gene}.log"
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

# Rule 6-2: MAFFT Alignment
rule mafft_alignment:
    input:
        merged_fasta = "results/06_alignment/{type}/merged/{gene}.fasta"
    params:
        ref_path = lambda w: os.path.join(config["paths"]["ref_alignment_dir"], f"{w.gene}.fasta"),
        mafft_opts = config["params"]["mafft"]["alignment_option"]
    output:
        aln = "results/06_alignment/{type}/aligned/{gene}.aln"
    threads: 4
    log: "logs/6-2/mafft_{type}_{gene}.log"
    shell:
        """
        if [ -s "{input.merged_fasta}" ]; then
            
            # Case 1: Reference exists -> use --add
            if [ -f "{params.ref_path}" ]; then
                echo "Reference found: {params.ref_path}" > {log}
                
                mafft {params.mafft_opts} --add {input.merged_fasta} \
                      --reorder {params.ref_path} > {output.aln} 2>> {log}
            
            # Case 2: No reference -> de novo
            else
                echo "No reference found. Running de novo." >> {log}

                mafft {params.mafft_opts} \
                      {input.merged_fasta} > {output.aln} 2>> {log}
            fi
            
        else
            touch {output.aln}
        fi
        """

# Rule 6-3: TrimAl
rule trimal_trimming:
    input:
        aln = "results/06_alignment/{type}/aligned/{gene}.aln"
    output:
        trimmed = "results/06_alignment/{type}/trimmed/{gene}.trimmed.aln"
    params:
        trimal_opts = config["params"]["trimal"]
    log: "logs/6-3/trimal_{type}_{gene}.log"
    shell:
        """
        if [ -s "{input.aln}" ]; then
            trimal -in {input.aln} -out {output.trimmed} {params.trimal_opts} > {log} 2>&1
        else
            touch {output.trimmed}
        fi
        """