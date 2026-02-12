# Phase 7: Phylogeny Tree Building (IQ-TREE)

# -----------------------------------------------------------------------------
# Rule: IQ-TREE
# Builds phylogenetic trees from trimmed alignments.
# Uses parameters defined in config/config.yaml (Model, Bootstrap, etc.)
# -----------------------------------------------------------------------------
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
    log: "logs/7/iqtree_{type}_{gene}.log"
    shell:
        """
        if [ -s "{input.aln}" ]; then
            
            # -s: Input alignment
            # --prefix: Output file prefix
            # -m: Model selection (MFP = ModelFinder Plus)
            # -B: Ultrafast Bootstrap iterations
            # -nt: Number of threads
            # -redo: Overwrite existing files
            
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