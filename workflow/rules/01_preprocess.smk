# workflow/rules/preprocessing.smk
# Rename header of fasta files and make blast database for each sample

RESULTS_DIR = config["paths"]["results"]    # path for output files (result directory)
GENOMES_DIR = config["paths"]["genomes"]    # path for input genome fasta files 

rule rename_fasta:
    input:
        fasta = lambda wildcards: f"{GENOMES_DIR}/{samples_df.loc[wildcards.sample_id, 'genome_filename']}"
    output:
        fasta = RESULTS_DIR + "/00_renamed/{sample_id}_renamed.fasta",
        map = RESULTS_DIR + "/00_renamed/{sample_id}_mapping.tsv"
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
        fna = RESULTS_DIR + "/00_renamed/{sample_id}_renamed.fasta"
    output:
        db_file = RESULTS_DIR + "/01_blastdb/{sample_id}_renamed.nhr"
    params:
        db_prefix = RESULTS_DIR + "/01_blastdb/{sample_id}_renamed"
    log:
        "logs/1-2/{sample_id}.log"
    shell:
        """
        makeblastdb -in {input.fna} -dbtype nucl -out {params.db_prefix} > {log} 2>&1
        """