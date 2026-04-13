#!/bin/bash
#SBATCH --job-name=aster_hypnales         # Job name
#SBATCH --mail-type=ALL               # Email notifications
#SBATCH --mail-user=seyeonkim@ufl.edu # Replace with your email
#SBATCH --nodes=1                     # Use 1 node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=2             # Number of CPU cores
#SBATCH --mem=4gb                     # Job memory
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/aster_hypnales%j.out     # Standard output and error log
#SBATCH --account=mcdaniel
#SBATCH --qos=mcdaniel-b

# Select "Order" you want to analyze (e.g., "Hypnales", "Pottiales", "Orthotrichales")
ORDER="Hypnales"

# 0. Set working directory (Modify if necessary)
## Input directories (iqtree results, (*.treefile))
IQTREE_RESULTS_DIR1=results/07_phylogeny/Best
IQTREE_RESULTS_DIR2=results/07_phylogeny/Not_sex_linked
## Metadata file: config/sample.tsv
METADATA_FILE=config/sample.tsv
## Output directory
OUTPUT_DIR=results/${ORDER}_aster/
### Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}
## Concatanated gene tree file name
CAT_TREE_FILE=results/${ORDER}/${ORDER}_concatenated.treefile

# 1. Creae taxa list for the selected order

# 2. Concatenate gene trees for the selected order

# Load modules
module load aster/1.22
# Run aster
wastral -i ${CAT_TREE_FILE} -o ${OUTPUT_DIR}/aster_output.tre -t 8


