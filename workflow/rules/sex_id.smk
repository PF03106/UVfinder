# BLAST for male markers
rule blast_male:
    input:
        query = config["paths"]["male_markers"],
        db_index = "results/01_blastdb/{sample_id}_renamed.nhr"
    output: "results/02_sex_id/{sample_id}_male.tblastn"
    params: db_prefix = "results/01_blastdb/{sample_id}_renamed", evalue = config["params"]["uv_blast_evalue"]
    log: "logs/2-1/blast_male_{sample_id}.log" 
    shell: "tblastn -query {input.query} -db {params.db_prefix} -evalue {params.evalue} -outfmt '6 std qlen' -out {output}"

# BLAST for female markers
rule blast_female:
    input:
        query = config["paths"]["female_markers"],
        db_index = "results/01_blastdb/{sample_id}_renamed.nhr"
    output: "results/02_sex_id/{sample_id}_female.tblastn"
    params: db_prefix = "results/01_blastdb/{sample_id}_renamed", evalue = config["params"]["uv_blast_evalue"]
    log: "logs/2-1/blast_female_{sample_id}.log"
    shell: "tblastn -query {input.query} -db {params.db_prefix} -evalue {params.evalue} -outfmt '6 std qlen' -out {output}"

# ID U or V or Unknown
rule assign_sex_differential:
    input:
        male_res = "results/02_sex_id/{sample_id}_male.tblastn",
        female_res = "results/02_sex_id/{sample_id}_female.tblastn"
    output: tsv = "results/02_sex_id/{sample_id}_sex_assignment.tsv"
    params:
        min_cov = config["params"]["min_coverage"],
        max_evalue = config["params"]["uv_blast_evalue"]
    log: "logs/2-2/assign_sex_{sample_id}.log"
    shell:
        """
        python3 workflow/scripts/assign_sex_diff.py \
            --male {input.male_res} \
            --female {input.female_res} \
            --output {output.tsv} \
            --min_cov {params.min_cov} \
            --max_evalue {params.max_evalue}
        """
# Aggregate all the result from sex_sgginment.tsv
rule aggregate_sex_id:
    input:
        results = expand("results/02_sex_id/{sample_id}_sex_assignment.tsv", sample_id=SAMPLES),
        metadata = "config/samples.tsv"
    output:
        all_res = "results/02_sex_id/all_samples_sex_assignment.tsv"
    log: "logs/2-3/aggregate_sex_id.log"
    run:
        import pandas as pd
        import os

        dfs=[]
        for f in input.results:
            df = pd.read_csv(f, sep='\t')
            sample_name = os.path.basename(f).replace("_sex_assignment.tsv", "")
            df.insert(0, "sample_id", sample_name)
            dfs.append(df)
        combined_df = pd.concat(dfs, ignore_index=True)

        meta_df = pd.read_csv(input.metadata, sep='\t')
        final_df = pd.merge(combined_df, meta_df[['sample_id', 'genus', 'species']], on='sample_id', how='left')
        cols = ['sample_id', 'genus', 'species'] + [col for col in final_df.columns if col not in ['sample_id', 'genus', 'species']]
        final_df = final_df[cols]
        final_df.to_csv(output.all_res, sep='\t', index=False)