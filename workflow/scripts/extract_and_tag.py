import pandas as pd
from Bio import SeqIO
import argparse
import os
import sys

def load_sex_linked_map(file_path):
    """
    sex_assign.tsv를 읽어 성염색체 번호와 성별을 매핑합니다.
    [수정됨] Unknown도 명시적으로 맵에 등록합니다.
    """
    linked_map = {}
    if not os.path.exists(file_path):
        return linked_map

    try:
        df = pd.read_csv(file_path, sep='\t')
        for _, row in df.iterrows():
            c_id = str(row['contig_id'])
            sex = str(row['sex'])
            
            if c_id.lower() == 'none':
                continue
                
            # 마지막 2자리 숫자 추출 (예: tig00000014 -> 14)
            chr_num = c_id[-2:]
            
            # [중요] 대소문자 구분 없이 처리
            sex_lower = sex.lower()
            
            if sex_lower == 'u':
                linked_map[chr_num] = 'U-linked'
            elif sex_lower == 'v':
                linked_map[chr_num] = 'V-linked'
            elif sex_lower == 'unknown': 
                linked_map[chr_num] = 'Unknown' # Unknown 상태 보존
            # 그 외(Autosome 등)는 맵에 없으면 나중에 Default 'A'로 처리됨
            
    except Exception as e:
        print(f"Warning parsing sex map: {e}")
            
    return linked_map

def extract_and_tag(blast_file, linked_map, genome_dict, out_dir, sample_id, mode='A'):
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen"]
    
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        return

    df = pd.read_csv(blast_file, sep='\t', names=cols)

    # [Start] 추가할 디버깅 코드 -------------------------
    print("\n--- 🔍 DEBUGGING START ---")
    print(f"1. BLAST에서 본 컨티그 이름 (상위 3개): {df['sseqid'].head(3).tolist()}")
    print(f"2. 게놈 파일의 컨티그 이름 (상위 3개): {list(genome_dict.keys())[:3]}")
    
    # 첫 번째 BLAST 결과가 게놈 사전에 있는지 테스트
    first_hit = str(df.iloc[0]['sseqid'])
    if first_hit in genome_dict:
        print(f"✅ 매칭 성공! '{first_hit}'을 찾았습니다.")
    else:
        print(f"❌ 매칭 실패! '{first_hit}'이 게놈 사전에 없습니다.")
        print("👉 띄어쓰기, 파이프(|), 버전 번호(.1) 등을 확인해보세요.")
    print("--- DEBUGGING END ---\n")
    # [End] ---------------------------------------------
    
    # B 모드(Best Hit)일 경우: 유전자(Locus)별로 가장 높은 점수 1개만 남김
    if mode == 'B':
        # G4471|... 형태에서 G4471만 임시 추출하여 중복 제거 기준(subset)으로 사용
        df['temp_locus_id'] = df['qseqid'].apply(lambda x: x.split('|')[0] if '|' in str(x) else x)
        df = df.sort_values('bitscore', ascending=False).drop_duplicates('temp_locus_id')

    for _, row in df.iterrows():
        raw_qseqid = str(row['qseqid'])
        sseqid = str(row['sseqid'])
        
        # [핵심 1] Qseqid 파싱 (G4471|OriginalID -> G4471 추출)
        if '|' in raw_qseqid:
            locus_id = raw_qseqid.split('|')[0] # G4471
        else:
            locus_id = raw_qseqid # 포맷이 다르면 그대로 사용

        # [핵심 2] 성염색체 태깅 (Unknown 처리 포함)
        hit_chr_num = sseqid[-2:]
        sex_type = linked_map.get(hit_chr_num, 'A') # 맵에 없으면 Autosome
        
        if sex_type == 'U-linked':
            tag_suffix = "U"
        elif sex_type == 'V-linked':
            tag_suffix = "V"
        elif sex_type == 'Unknown':
            tag_suffix = "Unk" # Unknown은 Unk 태그 부착
        else:
            tag_suffix = "A"
        
        # [핵심 3] 최종 헤더 구성: >Sample_Locus_Tag (통합 분석용)
        # 예: >S0035_G4471_U
        new_header_id = f"{sample_id}_{locus_id}_{tag_suffix}"
        
        # 서열 추출 (역상보 보정)
        start, end = int(row['sstart']), int(row['send'])
        is_reverse = start > end
        actual_start, actual_end = (end - 1, start) if is_reverse else (start - 1, end)
        
        if sseqid in genome_dict:
            seq_record = genome_dict[sseqid][actual_start:actual_end]
            if is_reverse:
                seq_record.seq = seq_record.seq.reverse_complement()
            
            # 레코드 정보 업데이트
            seq_record.id = new_header_id
            seq_record.description = f"original:{sseqid} qseqid:{raw_qseqid}"
            
            # [핵심 4] 파일 저장: 유전자 이름(G4471.fasta)으로 저장
            locus_file_path = os.path.join(out_dir, f"{locus_id}.fasta")
            
            with open(locus_file_path, "a") as f:
                SeqIO.write(seq_record, f, "fasta")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True, help="Sample ID (e.g., S0035)")
    parser.add_argument("--blast", required=True)
    parser.add_argument("--sex_map", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--out_best", required=True) # B_Best
    parser.add_argument("--out_all", required=True)  # A_all
    args = parser.parse_args()

    os.makedirs(args.out_best, exist_ok=True)
    os.makedirs(args.out_all, exist_ok=True)

    # 1. 매핑 정보 로드
    linked_map = load_sex_linked_map(args.sex_map)
    
    # 2. 유전체 로드
    genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

    # 3. 추출 및 태깅 실행
    extract_and_tag(args.blast, linked_map, genome_dict, args.out_all, args.sample, mode='A')
    extract_and_tag(args.blast, linked_map, genome_dict, args.out_best, args.sample, mode='B')

if __name__ == "__main__":
    main()