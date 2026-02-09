# Phase 6: Subsetting, Trimming, and Tree Inference

rule subset_samples_only:
    """
    Step 6.1: Remove reference sequences from the alignment.
    Keep only samples listed in samples.tsv.
    """
    input:
        aligned = "results/05_alignment/final_msa/Global_{locus}_aligned.fasta",
        samples_tsv = "config/samples.tsv"
    output:
        subset = "results/06_phylogeny/subset/Global_{locus}_samples.fasta"
    shell:
        """
        python3 workflow/scripts/subset_alignment.py \
            --input {input.aligned} \
            --output {output.subset} \
            --samples_tsv {input.samples_tsv}
        """

rule trimming_with_trimal:
    """
    Step 6.2: Trim the subset alignment (Samples Only).
    """
    input:
        # [수정] 입력이 전체 정렬이 아니라 'subset' 결과입니다.
        aligned = rules.subset_samples_only.output.subset
    output:
        trimmed = "results/06_phylogeny/trimmed/Global_{locus}_trimmed.fasta",
        html_log = "results/06_phylogeny/trimmed/Global_{locus}_report.html"
    log: "logs/phase6_trimal_{locus}.log"
    shell:
        """
        # -automated1: 샘플들끼리의 정렬 상태를 보고 최적의 트리밍 수행
        trimal -in {input.aligned} \
               -out {output.trimmed} \
               -htmlout {output.html_log} \
               -automated1 \
               > {log} 2>&1
        """

rule iqtree_inference:
    """
    Step 6.3: Infer Tree using IQ-TREE (Samples Only).
    """
    input:
        trimmed = rules.trimming_with_trimal.output.trimmed
    output:
        tree = "results/06_phylogeny/trees/Global_{locus}.treefile",
        iqtree_log = "results/06_phylogeny/trees/Global_{locus}.log"
    params:
        out_prefix = "results/06_phylogeny/trees/Global_{locus}"
    threads: 8 # (또는 4, 상황에 맞게)
    log: "logs/phase6_iqtree_{locus}.log"
    shell:
        """
        iqtree -s {input.trimmed} \
               -st AA \
               -m TEST \
               -bb 1000 \
               -nt {threads} \
               -pre {params.out_prefix} \
               > {log} 2>&1
        """