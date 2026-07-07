# workflow/rules/02_assign_sex.smk
# tBLASTn sex marker genes against each sample's blast database and assign sex and sex chromosome number.

RESULTS_DIR = config["paths"]["results"]    # path for output files (result directory)
SAMPLES_PATH = config["paths"]["samples_tsv"]

# BLAST for male markers
rule blast_male:
    input:
        query = config["paths"]["male_markers"],
        db_index = f"{RESULTS_DIR}/01_blastdb/{{sample_id}}_renamed.nhr"
    output: f"{RESULTS_DIR}/02_sex_id/{{sample_id}}_male.tblastn"
    params: db_prefix = f"{RESULTS_DIR}/01_blastdb/{{sample_id}}_renamed", evalue = config["params"]["uv_blast_evalue"]
    log: f"logs/2-1/blast_male_{{sample_id}}.log"
    shell: 
        """
        tblastn \
        -query {input.query} \
        -db {params.db_prefix} \
        -evalue {params.evalue} \
        -outfmt '6 std' -out {output} > {log} 2>&1
        """

# BLAST for female markers
rule blast_female:
    input:
        query = config["paths"]["female_markers"],
        db_index = f"{RESULTS_DIR}/01_blastdb/{{sample_id}}_renamed.nhr"
    output: f"{RESULTS_DIR}/02_sex_id/{{sample_id}}_female.tblastn"
    params: db_prefix = f"{RESULTS_DIR}/01_blastdb/{{sample_id}}_renamed", evalue = config["params"]["uv_blast_evalue"]
    log: f"logs/2-1/blast_female_{{sample_id}}.log"
    shell: 
        """
        echo "[$(date)] Starting BLAST search for {wildcards.sample_id}..." > {log}
        tblastn \
        -query {input.query} \
        -db {params.db_prefix} \
        -evalue {params.evalue} \
        -outfmt '6 std' -out {output} \
        >> {log} 2>&1
        echo "[$(date)] Finished BLAST search for {wildcards.sample_id}." >> {log}
        """

# ID U or V or Unknown
rule assign_sex_differential:
    input:
        male_res = f"{RESULTS_DIR}/02_sex_id/{{sample_id}}_male.tblastn",
        female_res = f"{RESULTS_DIR}/02_sex_id/{{sample_id}}_female.tblastn",
        metadata = f"{SAMPLES_PATH}",
        m_marker = config["paths"]["male_markers"], 
        f_marker = config["paths"]["female_markers"],  
    output: out_file = f"{RESULTS_DIR}/02_sex_id/{{sample_id}}_sex_assignment.tsv"
    params:
        min_bitscore_ratio_UV = config["params"]["min_bitscore_ratio_UV"]
    log: f"logs/2-2/assign_sex_{{sample_id}}.log"
    shell:
        """
        echo "[$(date)] working on sample {wildcards.sample_id}..." > {log}
        python3 workflow/scripts/assign_sex_diff.py \
            --sample {wildcards.sample_id} \
            --samples_tsv {input.metadata} \
            --male_blast {input.male_res} \
            --female_blast {input.female_res} \
            --male_marker {input.m_marker} \
            --female_marker {input.f_marker} \
            --output {output.out_file} \
            --min_bitscore_ratio_UV {params.min_bitscore_ratio_UV} >> {log} 2>&1
        """

# Aggregate all the result from sex_sgginment.tsv
rule aggregate_sex_id:
    input:
        results = expand(f"{RESULTS_DIR}/02_sex_id/{{sample_id}}_sex_assignment.tsv", sample_id=SAMPLES),
        metadata = f"{SAMPLES_PATH}",
    output:
        all_res = f"{RESULTS_DIR}/02_sex_id/all_samples_sex_assignment.tsv"
    log: f"logs/2-3/aggregate_sex_id.log"
    shell: 
        """
        echo "[$(date)] Run workflow/scripts/aggregate_sex_id.py"
        python3 -u workflow/scripts/aggregate_sex_id.py \
            --inputs {input.results} \
            --metadata {input.metadata} \
            --output {output.all_res} > {log} 2>&1
        echo "[$(date)] Finished aggregating tsv. Check {output.all_res} for results."
        """