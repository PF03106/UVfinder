# Phase 7: Phylogeny Tree Building (IQ-TREE)
# Rule 7.1: Build Tree with IQ-TREE
rule iqtree_build:
    input:
        aln = "results/06_alignment/{type}/trimmed/{gene}.trimmed.aln"
    output:
        treefile = "results/07_phylogeny/{type}/{gene}.treefile",
        report   = "results/07_phylogeny/{type}/{gene}.iqtree"
    params:
        prefix = "results/07_phylogeny/{type}/{gene}",
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
        treefile = "results/07_phylogeny/Best/{gene}.treefile", # <-- 여기에 콤마 추가!
        metadata = "config/samples.tsv"
    output:
        viz = "results/07_phylogeny/Best/{gene}_viz.png"
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
