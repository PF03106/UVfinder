#!/usr/bin/env python3
import os
import re
import sys

def clean_original_id(header_line):
    """
    1. Chromosome/LG 패턴이 명확하면 짧게 줄임 (예: Chr1, ChrX)
    2. 그 외(Scaffold 등)는 원래 헤더 전체를 가져와서 특수문자만 '_'로 교체함 (안전 제일)
    """
    # '>' 제거 및 앞뒤 공백 제거
    full_header = header_line.strip()[1:]
    
    # --- [1] Chromosome / Linkage Group 패턴 처리 (이건 줄이는 게 이득) ---
    # 예: "Chromosome 1", "chrX", "LG02" 등 숫자나 외자 알파벳(X,Y)이 뒤에 오는 경우만 잡음
    chrom_match = re.search(r'(?:chromosome|chr|LG)\s*[-_]?\s*([0-9]+|[XYZW]|Un)\b', full_header, re.IGNORECASE)
    
    if chrom_match:
        # Chr1, ChrX, ChrUn 등으로 변환
        return f"Chr{chrom_match.group(1)}"

    # --- [2] 그 외 모든 경우 (Scaffold 포함) : 안전하게 전체 보존 ---
    # 요청하신 대로 긴 이름을 그대로 살리되, 파일 시스템에 문제없게만 처리합니다.
    
    # 1. 공백, 탭, 쉼표, 파이프(|)를 모두 언더바(_)로 변경
    cleaned = re.sub(r'[ \t,;|]+', '_', full_header)
    
    # 2. 알파벳, 숫자, 점(.), 언더바(_), 하이픈(-)을 제외한 특수문자 제거 (괄호 등)
    cleaned = re.sub(r'[^a-zA-Z0-9._-]', '', cleaned)
    
    # 3. 언더바가 연속으로 나오면 하나로 합침 (Scaffold__14 -> Scaffold_14)
    cleaned = re.sub(r'_+', '_', cleaned).strip('_')
    
    return cleaned

def rename_fasta_headers(input_fna, output_fna, mapping_fna, s_id, order, genus, species):
    # 파이프라인 표준 접두사 (이미 길기 때문에 뒤에 붙는 ID가 중요함)
    prefix = f"{s_id}_{order}_{genus}_{species}"
    
    with open(input_fna, "r", encoding="utf-8", errors="ignore") as f_in, \
         open(output_fna, "w", encoding="utf-8") as f_out, \
         open(mapping_fna, "w", encoding="utf-8") as f_map:
        
        f_map.write("original_id\tnew_id\n")
        
        for line in f_in:
            if line.startswith(">"):
                original_header = line.strip()
                
                # 수정된 로직 적용
                new_clean_id = clean_original_id(original_header)
                
                # 최종 헤더: >S0029_..._MU104855.1_Ceratodon_purpureus_..._scaffold_14
                new_header = f">{prefix}_{new_clean_id}"
                
                f_out.write(new_header + "\n")
                # 매핑 로그 저장
                f_map.write(f"{original_header[1:]}\t{new_header[1:]}\n")
            else:
                f_out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python3 rename_fasta_headers.py <in.fna> <out.fna> <map.tsv> <ID> <Order> <Genus> <Species>")
        sys.exit(1)

    input_path, output_path, mapping_path = sys.argv[1:4]
    s_id, order_name, genus_name, species_name = sys.argv[4:8]

    print(f"🚀 Processing: {s_id} - Safe renaming mode (Keep full headers)")
    rename_fasta_headers(input_path, output_path, mapping_path, s_id, order_name, genus_name, species_name)
    print(f"✅ Renaming complete. Log saved to {mapping_path}")