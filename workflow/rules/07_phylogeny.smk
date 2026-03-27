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
    threads: config["params"]["iqtree"]["threads"]
    log: "logs/7-1/iqtree_{type}_{gene}.log"
    shell:
        """
        if [ -s "{input.aln}" ]; then
            iqtree -s {input.aln} \
                   --prefix {params.prefix} \
                   -m {params.model} \
                   -B {params.bootstrap} \
                   -nt {threads} \
                   -redo \
                   > {log} 2>&1
        else
            echo "Empty alignment file. Skipping Tree building." > {log}
            touch {output.treefile} {output.report}
        fi
        """

# Rule 7.2: Plot Tree with ETE Toolkit (Only for 'Best' hits)
rule plot_tree:
    input:
        treefile = f"{RESULTS_DIR}/07_phylogeny/Best/{{gene}}.treefile",
        metadata = "config/samples.tsv"
    output:
        viz = f"{RESULTS_DIR}/07_phylogeny/Best/{{gene}}_viz.png"
    log:
        "logs/7-2/plot_tree_Best_{gene}.log"
    shell:
        """
        if [ -s "{input.treefile}" ]; then
            python3 workflow/scripts/plot_tree.py \
                --tree {input.treefile} \
                --metadata {input.metadata} \
                --output_prefix {output.viz} \
                > {log} 2>&1
        else
            echo "Tree file not found. Skipping visualization." > {log}
            touch {output.viz}
        fi
        """
