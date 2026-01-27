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
rule make_blast_db:
    input:
        fna = "results/00_renamed/{sample_id}_renamed.fasta"
    output:
        db_file = "results/00_renamed/{sample_id}_renamed.fasta.nhr"
    params:
        db_name = "results/00_renamed/{sample_id}_renamed.fasta"
    log:
        "logs/blastdb/{sample_id}.log"
    shell:
        """
        makeblastdb -in {input.fna} -dbtype nucl -out {params.db_name} > {log} 2>&1
        """