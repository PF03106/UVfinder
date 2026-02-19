# Phase 7: Phylogeny Tree Building (IQ-TREE)

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

# Rule 7.2: Plot Tree with ETE Toolkit
rule plot_tree:
    input:
        treefile = "results/07_phylogeny/{type}/{gene}.treefile"
    output:
        viz = "results/07_phylogeny/{type}/{gene}_viz.png"
    log:
        "logs/7-2/plot_tree_{type}_{gene}.log"
    shell:
        """
        if [ -s "{input.treefile}" ]; then
            python3 workflow/scripts/plot_tree.py \
                --tree {input.treefile} \
                --output_prefix {output.viz} \
                > {log} 2>&1
        else
            echo "Tree file not found. Skipping visualization." > {log}
            touch {output.viz}
        fi
        """