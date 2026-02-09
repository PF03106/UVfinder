import os
import re
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def load_sex_linked_map(file_path):
    """
    성염색체 매핑 파일(TSV) 로드
    """
    linked_map = {}
    try:
        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            return linked_map

        df = pd.read_csv(file_path, sep='\t')
        for _, row in df.iterrows():
            c_id = str(row['contig_id'])
            sex = str(row['sex'])
            if c_id.lower() == 'none': continue
            
            match = re.search(r'[Cc]hr([0-9]+|[A-Za-z]+)', c_id)
            chr_num = match.group(1) if match else c_id
            
            sex_lower = sex.lower()
            if sex_lower == 'u': linked_map[chr_num] = 'U-linked'
            elif sex_lower == 'v': linked_map[chr_num] = 'V-linked'
            elif sex_lower == 'unknown': linked_map[chr_num] = 'Unknown'
    except Exception as e:
        print(f"Warning parsing sex map: {e}")
    return linked_map

def get_rank_letter(n):
    if n < 50: return f"R{n + 1}"
    else: return None

def extract_and_tag(blast_file, linked_map, genome_dict, out_best_dir, out_all_dir, sample_id, flank_bp=20):
    """
    BLAST 결과를 파싱하여 유전자별로 개별 파일을 생성하여 저장
    """
    # 1. [매우 중요] 출력 디렉토리(폴더) 생성
    # Snakemake는 이 폴더가 존재해야 작업 완료로 간주합니다.
    os.makedirs(out_best_dir, exist_ok=True)
    os.makedirs(out_all_dir, exist_ok=True)
    
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen"]
    
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        print(f"Warning: BLAST file not found or empty: {blast_file}")
        # 빈 폴더만 생성되고 종료됨 (에러 아님)
        return

    try:
        df = pd.read_csv(blast_file, sep='\t', names=cols)
        # 정렬
        df = df.sort_values(by=['qseqid', 'evalue', 'bitscore'], ascending=[True, True, False])
        
        # 2. Query ID(유전자) 별로 그룹화하여 처리
        # 여기서 qseqid는 'G4471|Original_Header' 형태
        for qseqid_raw, group in df.groupby("qseqid"):
            
            # 유전자 ID 추출 (파일명으로 사용)
            # 예: G4471|... -> G4471
            clean_locus = qseqid_raw.split('|')[0]
            
            gene_records = []
            
            # 해당 유전자의 히트들 처리
            for i, (_, row) in enumerate(group.iterrows()):
                rank = get_rank_letter(i)
                if rank is None: continue # 50등 밖은 버림
                
                chrom = row['sseqid']
                if chrom not in genome_dict: continue
                    
                full_seq = genome_dict[chrom].seq
                chrom_len = len(full_seq)
                sstart, send = int(row['sstart']), int(row['send'])
                
                # 좌표 및 strand 처리
                if sstart < send:
                    strand = 1
                    start_pos = max(0, sstart - 1 - flank_bp)
                    end_pos = min(chrom_len, send + flank_bp)
                    seq_extracted = full_seq[start_pos:end_pos]
                else:
                    strand = -1
                    start_pos = max(0, send - 1 - flank_bp)
                    end_pos = min(chrom_len, sstart + flank_bp)
                    seq_extracted = full_seq[start_pos:end_pos].reverse_complement()
                
                # Linkage info
                match = re.search(r'[Cc]hr([0-9]+|[A-Za-z]+)', chrom)
                chr_key = match.group(1) if match else chrom
                linkage_info = linked_map.get(chr_key, "Autosomal")
                
                # 헤더 포맷: >S0049_G4471_R1
                new_id = f"{sample_id}_{clean_locus}_{rank}"
                description = f"[{linkage_info}] original:{chrom}:{sstart}-{send}({strand})"
                
                record = SeqRecord(seq_extracted, id=new_id, description=description)
                gene_records.append(record)
            
            # 3. 파일 저장 (유전자별 개별 파일)
            if gene_records:
                # 파일명: G4471.fasta
                file_name = f"{clean_locus}.fasta"
                
                # (1) A_all 폴더에 저장 (모든 히트)
                out_path_all = os.path.join(out_all_dir, file_name)
                with open(out_path_all, "w") as f_all:
                    SeqIO.write(gene_records, f_all, "fasta")
                
                # (2) B_Best 폴더에 저장 (R1만 필터링)
                best_records = [rec for rec in gene_records if rec.id.endswith('_R1')]
                if best_records:
                    out_path_best = os.path.join(out_best_dir, file_name)
                    with open(out_path_best, "w") as f_best:
                        SeqIO.write(best_records, f_best, "fasta")

        print(f"Processed BLAST results into directories: {out_all_dir}, {out_best_dir}")

    except Exception as e:
        print(f"Error processing BLAST file: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--blast", required=True)
    parser.add_argument("--sex_map", required=True)
    parser.add_argument("--genome", required=True)
    # 인자 이름 변경: directory임을 명시
    parser.add_argument("--out_best_dir", required=True, help="Output Directory for Best hits")
    parser.add_argument("--out_all_dir", required=True, help="Output Directory for All hits")
    
    args = parser.parse_args()

    print(f"Loading genome: {args.genome}...")
    genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    
    sex_map = load_sex_linked_map(args.sex_map)
    
    extract_and_tag(
        blast_file=args.blast,
        linked_map=sex_map,
        genome_dict=genome_dict,
        out_best_dir=args.out_best_dir,
        out_all_dir=args.out_all_dir,
        sample_id=args.sample
    )