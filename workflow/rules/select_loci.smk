rule select_interesting_loci:
    input:
        fasta_files = expand("results/03_locus_search/{sample}/A_all/{sample}_extracted.fasta", sample=SAMPLES)
    output:
        locus_list = "results/04_filtering/interesting_loci.txt"
    log:
        "logs/select_loci.log"
    script:
        "scripts/select_loci.py"