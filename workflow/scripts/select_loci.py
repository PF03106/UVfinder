import os
import argparse
from Bio import SeqIO

def select_all_sex_linked_loci(input_dirs, output_file):
    """
    Scans all provided directories and lists every locus that has
    at least one sequence tagged with _U or _V.
    """
    interesting_loci = set()
    
    for sample_dir in input_dirs:
        if not os.path.exists(sample_dir):
            continue
            
        # B_Best 폴더 내의 각 유전자 파일을 전수 조사합니다
        for fasta_name in os.listdir(sample_dir):
            if not fasta_name.endswith(".fasta"):
                continue
            
            locus_id = fasta_name.replace(".fasta", "")
            file_path = os.path.join(sample_dir, fasta_name)
            
            try:
                for record in SeqIO.parse(file_path, "fasta"):
                    # 헤더에서 정확히 _U 또는 _V가 포함되어 있는지 확인합니다
                    if "_U" in record.id or "_V" in record.id:
                        interesting_loci.add(locus_id)
                        break # 이 로커스는 이미 대상에 포함되었으므로 다음 파일로 이동
            except Exception as e:
                print(f"Error parsing {file_path}: {e}")

    # 추출된 로커스 ID를 알파벳 순으로 정렬하여 저장합니다
    with open(output_file, "w") as f:
        for locus in sorted(list(interesting_loci)):
            f.write(f"{locus}\n")
            
    print(f"✅ Census complete. {len(interesting_loci)} unique loci identified with sex-linked signals.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Census of all sex-linked loci across samples.")
    parser.add_argument("--input_dirs", nargs="+", required=True, help="Directories of B_Best hits")
    parser.add_argument("--output", required=True, help="Output list of locus IDs")
    
    args = parser.parse_args()
    select_all_sex_linked_loci(args.input_dirs, args.output)