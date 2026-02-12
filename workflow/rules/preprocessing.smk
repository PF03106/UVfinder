# workflow/rules/preprocessing.smk

rule rename_fasta:
    input:
        fasta = lambda wildcards: f"resources/genomes/{samples_df.loc[wildcards.sample_id, 'genome_filename']}"
    output:
        fasta = "results/00_renamed/{sample_id}_renamed.fasta",
        map = "results/00_renamed/{sample_id}_mapping.tsv"
    params:
        #  wildcards.sample_id from Snakefile SAMPLES list
        order   = lambda wildcards: samples_df.loc[wildcards.sample_id, "order"],
        genus   = lambda wildcards: samples_df.loc[wildcards.sample_id, "genus"],
        species = lambda wildcards: samples_df.loc[wildcards.sample_id, "species"]
    log:
        "logs/1-1/{sample_id}.log"
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
rule make_blast_db:
    input:
        fna = "results/00_renamed/{sample_id}_renamed.fasta"
    output:
        db_file = "results/01_blastdb/{sample_id}_renamed.nhr"
    params:
        db_prefix = "results/01_blastdb/{sample_id}_renamed"
    log:
        "logs/1-2/{sample_id}.log"
    shell:
        """
        mkdir -p results/01_blastdb
        makeblastdb -in {input.fna} -dbtype nucl -out {params.db_prefix} > {log} 2>&1
        """