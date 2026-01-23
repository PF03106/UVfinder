# workflow/rules/preprocessing.smk

rule rename_fasta:
    input:
        fasta = "resources/genomes/{sample_id}.fasta"
    output:
        fasta = "results/00_renamed/{sample_id}_renamed.fasta",
        map = "results/00_renamed/{sample_id}_mapping.tsv"
    params:
        #  wildcards.sample_id from Snakefile SAMPLES list
        order   = lambda wildcards: samples_df.loc[wildcards.sample_id, "order"],
        genus   = lambda wildcards: samples_df.loc[wildcards.sample_id, "genus"],
        species = lambda wildcards: samples_df.loc[wildcards.sample_id, "species"]
    log:
        "logs/rename/{sample_id}.log"
    shell:
        """
        python3 workflow/scripts/rename_fasta_headers.py \
            {input.fasta} \
            {output.fasta} \
            {output.map} \
            {wildcards.sample_id} \
            {params.order} \
            {params.genus} \
            {params.species} > {log} 2>&1
        """