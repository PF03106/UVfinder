import pandas as pd
from Bio import SeqIO
import argparse
import os
import sys

def load_sex_linked_map(file_path):
    """
    sex_assign.tsv를 읽어 성염색체 번호(마지막 2자리)와 성별을 매핑
    """
    linked_map = {}
    if not os.path.exists(file_path):
        return linked_map

    try:
        df = pd.read_csv(file_path, sep='\t')
        for _, row in df.iterrows():
            c_id = str(row['contig_id'])
            sex = str(row['sex'])
            
            if c_id.lower() == 'none' or sex.lower() == 'unknown':
                continue
                
            # 마지막 2자리 숫자 추출 (예: tig00000014 -> 14)
            chr_num = c_id[-2:]
            if sex == 'U':
                linked_map[chr_num] = 'U-linked'
            elif sex == 'V':
                linked_map[chr_num] = 'V-linked'
    except Exception as e:
        print(f"Warning parsing sex map: {e}")
            
    return linked_map

def extract_and_tag(blast_file, linked_map, genome_dict, out_dir, sample_id, mode='A'):
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen"]
    
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        return

    df = pd.read_csv(blast_file, sep='\t', names=cols)
    
    # B 모드(Best Hit)일 경우 Bit-score 기준 유전자별 Top hit만 남김
    if mode == 'B':
        # G4471|... 형태에서 G4471만 보고 중복 제거를 해야 함
        df['locus_id'] = df['qseqid'].apply(lambda x: x.split('|')[0] if '|' in str(x) else x)
        df = df.sort_values('bitscore', ascending=False).drop_duplicates('locus_id')

    for _, row in df.iterrows():
        raw_qseqid = str(row['qseqid'])
        sseqid = str(row['sseqid'])
        
        # [핵심 수정 1] Qseqid 파싱 (G4471|OriginalID -> G4471 추출)
        if '|' in raw_qseqid:
            locus_id = raw_qseqid.split('|')[0] # G4471
        else:
            locus_id = raw_qseqid # 만약 포맷이 안맞으면 그대로 사용

        # 1. 성염색체 태깅
        hit_chr_num = sseqid[-2:]
        sex_type = linked_map.get(hit_chr_num, 'A') # 기본값 A (Autosome)
        
        # 태그 포맷 결정 (_U, _V, _A)
        if sex_type == 'U-linked':
            tag_suffix = "U"
        elif sex_type == 'V-linked':
            tag_suffix = "V"
        else:
            tag_suffix = "A"
        
        # [핵심 수정 2] 최종 헤더 구성: >Sample_Locus_Tag
        # 예: >S0035_G4471_U
        new_header_id = f"{sample_id}_{locus_id}_{tag_suffix}"
        
        # 2. 서열 추출
        start, end = int(row['sstart']), int(row['send'])
        is_reverse = start > end
        actual_start, actual_end = (end - 1, start) if is_reverse else (start - 1, end)
        
        if sseqid in genome_dict:
            seq_record = genome_dict[sseqid][actual_start:actual_end]
            if is_reverse:
                seq_record.seq = seq_record.seq.reverse_complement()
            
            # 레코드 정보 업데이트
            seq_record.id = new_header_id
            seq_record.description = f"original_contig:{sseqid} qseqid:{raw_qseqid}"
            
            # 3. 파일 저장 (유전자별 파일 생성: G4471.fasta)
            # 폴더 안에 G4471.fasta가 생성되고 그 안에 해당 서열이 추가됨
            locus_file_path = os.path.join(out_dir, f"{locus_id}.fasta")
            
            with open(locus_file_path, "a") as f:
                SeqIO.write(seq_record, f, "fasta")

def main():
    parser = argparse.ArgumentParser()
    # [수정 3] --sample 인자 추가
    parser.add_argument("--sample", required=True, help="Sample ID (e.g., S0035)")
    parser.add_argument("--blast", required=True)
    parser.add_argument("--sex_map", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--out_best", required=True) # B_Best
    parser.add_argument("--out_all", required=True)  # A_all
    args = parser.parse_args()

    os.makedirs(args.out_best, exist_ok=True)
    os.makedirs(args.out_all, exist_ok=True)

    linked_map = load_sex_linked_map(args.sex_map)
    
    # 유전체 로드 (대용량일 경우 인덱싱 사용)
    genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

    # A와 B 실행 시 sample ID 전달
    extract_and_tag(args.blast, linked_map, genome_dict, args.out_all, args.sample, mode='A')
    extract_and_tag(args.blast, linked_map, genome_dict, args.out_best, args.sample, mode='B')

if __name__ == "__main__":
    main()