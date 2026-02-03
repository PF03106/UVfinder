# Male 마커 전용 BLAST
rule blast_male:
    input:
        query = config["paths"]["male_markers"],
        db_index = "results/01_blastdb/{sample_id}_renamed.nhr"
    output: "results/02_sex_id/{sample_id}_male.tblastn"
    params: db_prefix = "results/01_blastdb/{sample_id}_renamed", evalue = config["params"]["uv_blast_evalue"]
    shell: "tblastn -query {input.query} -db {params.db_prefix} -evalue {params.evalue} -outfmt '6 std qlen' -out {output}"

# Female 마커 전용 BLAST
rule blast_female:
    input:
        query = config["paths"]["female_markers"],
        db_index = "results/01_blastdb/{sample_id}_renamed.nhr"
    output: "results/02_sex_id/{sample_id}_female.tblastn"
    params: db_prefix = "results/01_blastdb/{sample_id}_renamed", evalue = config["params"]["uv_blast_evalue"]
    shell: "tblastn -query {input.query} -db {params.db_prefix} -evalue {params.evalue} -outfmt '6 std qlen' -out {output}"

# 두 결과를 비교하여 성별 판정
rule assign_sex_differential:
    input:
        male_res = "results/02_sex_id/{sample_id}_male.tblastn",
        female_res = "results/02_sex_id/{sample_id}_female.tblastn"
    output: tsv = "results/02_sex_id/{sample_id}_sex_assignment.tsv"
    params:
        min_cov = config["params"]["min_coverage"],
        min_id = config["params"]["min_identity"]
    shell:
        """
        python3 workflow/scripts/assign_sex_diff.py \
            --male {input.male_res} \
            --female {input.female_res} \
            --output {output.tsv} \
            --min_cov {params.min_cov} \
            --min_id {params.min_id}
        """