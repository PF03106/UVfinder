import pandas as pd
from Bio import SeqIO
import argparse
import os

def load_sex_linked_map(file_path):
    """
    sex_assign.tsv를 읽어 성염색체 번호(마지막 2자리)와 성별을 매핑합니다.
    예: {'03': 'U-linked', '14': 'V-linked'}
    """
    linked_map = {}
    if not os.path.exists(file_path):
        return linked_map

    df = pd.read_csv(file_path, sep='\t')
    
    for _, row in df.iterrows():
        c_id = str(row['contig_id'])
        sex = str(row['sex'])
        
        # 'none'이거나 'Unknown'이면 스킵
        if c_id.lower() == 'none' or sex.lower() == 'unknown':
            continue
            
        # 마지막 2자리 숫자 추출
        chr_num = c_id[-2:]
        if sex == 'U':
            linked_map[chr_num] = 'U-linked'
        elif sex == 'V':
            linked_map[chr_num] = 'V-linked'
            
    return linked_map

def extract_and_tag(blast_file, linked_map, genome_dict, out_dir, mode='A'):
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen"]
    
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        return

    df = pd.read_csv(blast_file, sep='\t', names=cols)
    
    # B 모드(Best Hit)일 경우 Bit-score 기준 Top hit만 남김
    if mode == 'B':
        df = df.sort_values('bitscore', ascending=False).drop_duplicates('qseqid')

    for _, row in df.iterrows():
        locus = row['qseqid']
        sseqid = row['sseqid']
        
        # 1. 성염색체 연관 태깅 결정 (마지막 2자리 비교)
        hit_chr_num = sseqid[-2:]
        tag_suffix = linked_map.get(hit_chr_num, 'A')
        tag = f"_{tag_suffix}"
        
        # 2. 서열 추출 (역상보 보정 포함)
        start, end = int(row['sstart']), int(row['send'])
        is_reverse = start > end
        actual_start, actual_end = (end - 1, start) if is_reverse else (start - 1, end)
        
        if sseqid in genome_dict:
            seq_record = genome_dict[sseqid][actual_start:actual_end]
            if is_reverse:
                seq_record.seq = seq_record.seq.reverse_complement()
            
            # 3. 헤더 형식: OriginalID_Tag
            seq_record.id = f"{sseqid}{tag}"
            seq_record.description = f"locus:{locus}"
            
            # 4. 파일 저장
            locus_file = os.path.join(out_dir, f"{locus}.fasta")
            with open(locus_file, "a") as f:
                SeqIO.write(seq_record, f, "fasta")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast", required=True)
    parser.add_argument("--sex_map", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--out_best", required=True) # B_Best
    parser.add_argument("--out_all", required=True)  # A_all
    args = parser.parse_args()

    os.makedirs(args.out_best, exist_ok=True)
    os.makedirs(args.out_all, exist_ok=True)

    # 1. 성염색체 번호 매핑 데이터 로드
    linked_map = load_sex_linked_map(args.sex_map)
    
    # 2. 유전체 인덱싱 (메모리 효율적 관리)
    genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

    # 3. A (All hits) 및 B (Best hits) 추출
    extract_and_tag(args.blast, linked_map, genome_dict, args.out_all, mode='A')
    extract_and_tag(args.blast, linked_map, genome_dict, args.out_best, mode='B')

if __name__ == "__main__":
    main()