# Phase 5: Global Alignment (All samples + All references)

rule order_based_merging:
    """
    Step 5.1: Merge ALL samples and ALL references into one global file per locus.
    Output: results/05_alignment/merged_with_ref/Global_{Locus}.fasta
    """
    input:
        loci_list = "results/04_filtered/interesting_loci.txt",
        sample_sequences = "results/04_filtered/collected_loci",
        ref_probes = "resources/query_sets" # 실제 폴더명 확인 완료
    output:
        merged_dir = directory("results/05_alignment/merged_with_ref")
    log: "logs/phase5_merging.log"
    shell:
        """
        python3 workflow/scripts/merge_with_ref.py \
            --loci_list {input.loci_list} \
            --sample_dir {input.sample_sequences} \
            --ref_dir {input.ref_probes} \
            --output_dir {output.merged_dir}
        """

rule mafft_add_alignment:
    """
    Step 5.2: Align the Global file using MAFFT.
    """
    input:
        # 파일명 패턴 변경: Global_{locus}.fasta
        merged = "results/05_alignment/merged_with_ref/Global_{locus}.fasta"
    output:
        aligned = "results/05_alignment/final_msa/Global_{locus}_aligned.fasta"
    log: "logs/phase5_mafft_{locus}.log"
    threads: 8 # 데이터가 커졌으니 쓰레드를 늘리는 게 좋습니다.
    shell:
        """
        # --auto 옵션을 사용하여 데이터 크기에 맞춰 알고리즘 자동 선택
        mafft --amino --auto --thread {threads} {input.merged} > {output.aligned} 2> {log}
        """