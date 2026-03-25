# workflow/rules/goflag_blast_locus_search.smk
# BLAST goflag probes against each sample's genome and devide the results into "Best" and "All" hits

RESULTS_DIR = config["paths"]["results"]    # path for output files (result directory)

# 1. Call Python script to generate order-specific probes
rule prepare_order_probes:
    input:
        probe_dir = "resources/query_sets", 
        samples_tsv = "config/samples.tsv"
    output:
        #specific_probes = temp("resources/probes_by_sample/{sample_id}_specific_probes.fasta")     
        specific_probes = "resources/probes_by_sample/{sample_id}_specific_probes.fasta"
    log: "logs/3-1/prepare_probes_{sample_id}.log"
    shell:
        """
        python3 workflow/scripts/blast_by_order.py \
            --sample {wildcards.sample_id} \
            --samples_tsv {input.samples_tsv} \
            --probe_dir {input.probe_dir} \
            --output {output.specific_probes} \
            > {log} 2>&1
        """

# 2. Run BLAST with custom probes
rule blast_search:
    input:
        genome = f"{RESULTS_DIR}/00_renamed/{{sample_id}}_renamed.fasta",
        query = "resources/probes_by_sample/{sample_id}_specific_probes.fasta"
    output:
        blast_out = f"{RESULTS_DIR}/03_locus_search/{{sample_id}}/{{sample_id}}_blast_results.txt"
    params:
        evalue = config["params"]["query_blast_evalue"],
        outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    log: f"logs/3-2/blast_{{sample_id}}.log"
    shell:
        """
        tblastn \
            -query {input.query} \
            -subject {input.genome} \
            -evalue {params.evalue} \
            -outfmt "{params.outfmt}" \
            -out {output.blast_out} \
            2> {log}
        """

# 3. Partition BLAST results into All hits and Best hits. 
rule partition_blast_results:
    input:
        blast = f"{RESULTS_DIR}/03_locus_search/{{sample_id}}/{{sample_id}}_blast_results.txt",
        sex_map = f"{RESULTS_DIR}/02_sex_id/{{sample_id}}_sex_assignment.tsv"
        overlap_threshold = config["params"].get("overlap_threshold", 0.5)
    output:
        best_tsv = f"{RESULTS_DIR}/03_locus_search/{{sample_id}}/B_Best_hits.tsv",
        all_tsv = f"{RESULTS_DIR}/03_locus_search/{{sample_id}}/A_All_hits.tsv"
    log: f"logs/3-3/All_or_Best_partition_{{sample_id}}.log"
    shell:
        """
        python3 workflow/scripts/extract_best_all.py \
            --blast {input.blast} \
            --sex_map {input.sex_map} \
            --out_best {output.best_tsv} \
            --out_all {output.all_tsv} \
            --overlap_threshold {input.overlap_threshold} \
            > {log} 2>&1
        """