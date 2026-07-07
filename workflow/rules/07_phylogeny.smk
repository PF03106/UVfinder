# workflow/rules/07_phylogeny.smk
# Phylogeny Tree Building (IQ-TREE)

RESULTS_DIR = config["paths"]["results"]    # path for output files (result directory)

# Rule 7.1: Build Tree with IQ-TREE
rule iqtree_build:
    input:
        aln = f"{RESULTS_DIR}/06_alignment/{{type}}/trimmed/{{gene}}.trimmed.aln"
    output:
        treefile = f"{RESULTS_DIR}/07_phylogeny/{{type}}/{{gene}}.treefile",
        report   = f"{RESULTS_DIR}/07_phylogeny/{{type}}/{{gene}}.iqtree"
    params:
        prefix = f"{RESULTS_DIR}/07_phylogeny/{{type}}/{{gene}}",
        model = config["params"]["iqtree"]["model"],       
        bootstrap = config["params"]["iqtree"]["bootstrap"]
    log: "logs/7-1/iqtree_{type}_{gene}.log"
    shell:
        """
        if [ -s "{input.aln}" ]; then
            iqtree -s {input.aln} \
                   --prefix {params.prefix} \
                   -m {params.model} \
                   -B {params.bootstrap} \
                   -T AUTO \
                   > {log} 2>&1
        else
            echo "Empty alignment file. Skipping Tree building." > {log}
            touch {output.treefile} {output.report}
        fi
        """
