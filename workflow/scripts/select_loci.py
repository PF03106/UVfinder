import sys
import os

# Snakemake 객체를 통해 입력과 출력을 받아옵니다.
input_files = snakemake.input.fasta_files
output_file = snakemake.output.locus_list

interesting_loci = set()

try:
    for f in input_files:
        with open(f, 'r') as file:
            for line in file:
                if line.startswith(">"):
                    # 헤더 예시: >S0035_bryales_Anomobryum_..._G4471_U
                    parts = line.strip().split("_")
                    
                    # 설계하신 로직: 끝에서 두 번째가 Locus ID, 마지막이 성별 태그
                    locus_id = parts[-2]
                    sex_tag = parts[-1]

                    # U 또는 V 태그가 하나라도 있으면 선별
                    if sex_tag.upper() in ["U", "V"]:
                        interesting_loci.add(locus_id)

    # 결과 저장
    with open(output_file, 'w') as out:
        for locus in sorted(list(interesting_loci)):
            out.write(f"{locus}\n")

except Exception as e:
    # 에러 발생 시 로그에 기록 (선택 사항)
    with open(snakemake.log[0], "w") as log_file:
        log_file.write(str(e))
    sys.exit(1)