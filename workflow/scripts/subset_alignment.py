import argparse
from Bio import SeqIO
import pandas as pd

def subset_alignment(input_aln, output_aln, samples_tsv):
    # 1. 살려둘 샘플 ID 리스트 만들기
    df = pd.read_csv(samples_tsv, sep='\t')
    # 샘플 ID가 'S0002'라면, FASTA 헤더는 'S0002_...' 형식이므로 
    # 매칭을 위해 ID 리스트를 준비합니다.
    valid_ids = set(df['sample_id'].astype(str))

    # 2. 필터링 및 저장
    kept_sequences = []
    
    for record in SeqIO.parse(input_aln, "fasta"):
        # 헤더의 첫 부분(S0002) 추출
        sample_id = record.id.split("_")[0]
        
        # 내 샘플 리스트에 있는 경우만 유지
        if sample_id in valid_ids:
            kept_sequences.append(record)

    # 3. 결과 저장 (만약 샘플이 하나도 없으면 빈 파일 방지용 체크 가능)
    if kept_sequences:
        SeqIO.write(kept_sequences, output_aln, "fasta")
        print(f"✅ Subset complete: Kept {len(kept_sequences)} sequences.")
    else:
        print(f"⚠️ Warning: No matching samples found in {input_aln}")
        # 빈 파일이라도 생성해야 Snakemake 에러 방지
        open(output_aln, 'w').close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--samples_tsv", required=True)
    args = parser.parse_args()
    
    subset_alignment(args.input, args.output, args.samples_tsv)