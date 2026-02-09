import os
import argparse
import glob
from Bio import SeqIO

def get_best_translation(dna_seq):
    """ Select the best 3frame aa translation """
    best_aa = None
    min_stops = float('inf')
    best_len = 0
    
    for frame in range(3):
        seq_len = len(dna_seq)
        trim_len = (seq_len - frame) // 3 * 3
        if trim_len == 0: continue
        
        translated = dna_seq[frame:frame+trim_len].translate()
        stop_count = translated.count("*")
        if translated.endswith("*"): stop_count -= 1
            
        if stop_count < min_stops:
            min_stops = stop_count
            best_aa = translated
            best_len = len(translated)
        elif stop_count == min_stops:
            if len(translated) > best_len:
                best_aa = translated
                best_len = len(translated)
    return best_aa

def merge_all_references(loci_list_path, sample_dir, ref_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    with open(loci_list_path, 'r') as f:
        target_loci = [line.strip() for line in f if line.strip()]

    print(f"--- 🌍 Creating Global Alignments for {len(target_loci)} loci ---")

    for locus in target_loci:
        # 1. 내 샘플 모두 가져오기 (13종 전체)
        sample_file = os.path.join(sample_dir, f"{locus}_all.fasta")
        if not os.path.exists(sample_file):
            continue

        merged_records = []
        
        # 샘플 DNA -> AA 번역 및 추가
        for record in SeqIO.parse(sample_file, "fasta"):
            aa_seq = get_best_translation(record.seq)
            if aa_seq:
                record.seq = aa_seq
                record.description += " [sample_translated]"
                merged_records.append(record)

        # 2. 레퍼런스 모두 찾아오기 (Order 상관없이 싹 다)
        # resource/query_sets 안에 있는 해당 유전자(G####)가 포함된 모든 파일을 찾습니다.
        # 예: query_sets/Bryales/G4471.fasta, query_sets/Ditrichales/G4471.fasta 등
        # 혹은 파일명에 locus ID가 포함된 모든 패턴 검색
        ref_files = glob.glob(os.path.join(ref_dir, "**", f"*{locus}*.fasta"), recursive=True)
        
        for rf in ref_files:
            try:
                ref_recs = list(SeqIO.parse(rf, "fasta"))
                # 레퍼런스 ID에 출처(파일명)를 표시해주면 나중에 트리 볼 때 편함
                for r in ref_recs:
                    r.description += f" [Ref_Source:{os.path.basename(rf)}]"
                merged_records.extend(ref_recs)
            except:
                print(f"⚠️ Warning: Could not parse reference file {rf}")

        # 3. 통합 파일 저장
        if merged_records:
            # 파일명 앞에 Global을 붙여 구분
            out_filename = f"Global_{locus}.fasta"
            out_path = os.path.join(output_dir, out_filename)
            
            # ID 중복 제거 (필수)
            unique_records = {rec.id: rec for rec in merged_records}.values()
            
            SeqIO.write(unique_records, out_path, "fasta")

    print(f"✅ Global merging complete. Files saved to {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--loci_list", required=True)
    parser.add_argument("--sample_dir", required=True)
    parser.add_argument("--ref_dir", required=True)
    # samples_tsv는 이제 필요 없습니다 (Order로 안 나누니까요)
    parser.add_argument("--output_dir", required=True)
    
    args = parser.parse_args()
    
    merge_all_references(
        args.loci_list, 
        args.sample_dir, 
        args.ref_dir, 
        args.output_dir
    )