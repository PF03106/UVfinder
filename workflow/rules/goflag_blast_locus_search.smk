# --- Phase 3: Whole-Genome Locus Search (GoFlag Probes) ---

rule prepare_sample_queries:
    """
    Step 3.1: Extract probe sequences matching the Order.
    [수정됨] 헤더에 Locus ID(G####)를 추가하여 BLAST 후에도 유전자 정보를 유지함.
    """
    input:
        probe_dir = "resources/query_sets" 
    output:
        temp_fasta = temp("results/03_locus_search/{sample_id}/matched_queries.fasta")
    run:
        import os
        from Bio import SeqIO
        
        # 1. Sample의 Order 가져오기
        try:
            sample_order = samples_df.loc[wildcards.sample_id, "order"]
        except KeyError:
            print(f"❌ Error: Sample ID '{wildcards.sample_id}' not found in samples.tsv.")
            raise

        order_norm = str(sample_order).lower().replace(" ", "").strip()
        print(f"--- 🧬 Extracting GoFlag probes for {wildcards.sample_id} (Target Order: {order_norm}) ---")
        
        count = 0
        with open(output.temp_fasta, "w") as out_f:
            # 400개의 유전자 파일을 순회
            for locus_file in os.listdir(input.probe_dir):
                if not locus_file.endswith(".fasta"):
                    continue
                
                # [핵심] 파일 이름에서 유전자 ID 추출 (예: G4471.AA.fasta -> G4471)
                locus_id = locus_file.split(".")[0] 
                
                locus_path = os.path.join(input.probe_dir, locus_file)
                for record in SeqIO.parse(locus_path, "fasta"):
                    if record.id.lower().startswith(f"{order_norm}_"):
                        # [핵심] 헤더 변경: >G4471|Original_Header
                        # 이렇게 해야 BLAST 결과의 qseqid에 G4471이 남습니다.
                        record.id = f"{locus_id}|{record.id}"
                        record.description = record.id 
                        SeqIO.write(record, out_f, "fasta")
                        count += 1
        
        if count == 0:
            print(f"⚠️ Warning: No matching probes found for order '{order_norm}' in {wildcards.sample_id}.")

rule blast_queries:
    """
    Step 3.2: Perform tBLASTn
    """
    input:
        query = rules.prepare_sample_queries.output.temp_fasta,
        db_index = "results/01_blastdb/{sample_id}_renamed.nhr"
    output:
        "results/03_locus_search/{sample_id}_goflag.tblastn"
    params:
        db_prefix = "results/01_blastdb/{sample_id}_renamed",
        evalue = config["params"]["query_blast_evalue"]
    shell:
        """
        tblastn -query {input.query} \
                -db {params.db_prefix} \
                -evalue {params.evalue} \
                -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
                -out {output}
        """

rule filter_and_tag_hits:
    """
    Step 3.3: Parse BLAST results, tag sex-linkage, and extract sequences.
    """
    input:
        blast_res = rules.blast_queries.output,
        sex_map = "results/02_sex_id/{sample_id}_sex_assignment.tsv", # 경로 확인 필요
        genome = "results/00_renamed/{sample_id}_renamed.fasta"
    output:
        all_hits = directory("results/03_locus_search/{sample_id}/A_all"),
        best_hits = directory("results/03_locus_search/{sample_id}/B_Best")
    shell:
        """
        # [핵심] python 스크립트 실행 및 --sample 인자 전달
        python3 workflow/scripts/extract_and_tag.py \
            --sample {wildcards.sample_id} \
            --blast {input.blast_res} \
            --sex_map {input.sex_map} \
            --genome {input.genome} \
            --out_best {output.best_hits} \
            --out_all {output.all_hits}
        """