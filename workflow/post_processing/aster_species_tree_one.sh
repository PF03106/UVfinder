#!/bin/bash
#SBATCH --job-name=aster_all         # Job name
#SBATCH --mail-type=ALL               # Email notifications
#SBATCH --mail-user=seyeonkim@ufl.edu # Replace with your email
#SBATCH --nodes=1                     # Use 1 node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=8             # Number of CPU cores
#SBATCH --mem=8gb                     # Job memory
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/aster_all%j.out     # Standard output and error log
#SBATCH --account=mcdaniel
#SBATCH --qos=mcdaniel-b

# 1. Set variables and working directory (Modify if necessary)
OUT_FILE_PREFIX=all_sex_linked
THREAD=8
## Input directories (iqtree results, (*.treefile))
IQTREE_RESULTS_DIR1=results/07_phylogeny_ALL/Best
## Output directory
OUTPUT_DIR=results/${OUT_FILE_PREFIX}_aster/
### Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}
## Concatanated gene tree file name
RAW_CAT_TREE_FILE=${OUTPUT_DIR}/${OUT_FILE_PREFIX}_raw_concatenated.treefile
CLEAN_CAT_TREE_FILE=${OUTPUT_DIR}/${OUT_FILE_PREFIX}_clean_concatenated.treefile

# 2. Concatenate gene trees for the selected order
cat ${IQTREE_RESULTS_DIR1}/*.treefile > ${RAW_CAT_TREE_FILE}  # Modify if necessary
if [ $? -ne 0 ]; then
    echo "Error concatenating gene trees. Please check the input directories and files."
    exit 1
fi
echo "Concatenated gene trees saved to ${RAW_CAT_TREE_FILE}"

# 3. clean the concatenated tree file by removing redundant names
sed -E 's/_G[0-9]+_R[0-9]+_[AUVN]//g' ${RAW_CAT_TREE_FILE} > ${CLEAN_CAT_TREE_FILE}

# 4. Run WASTRAL(Aster)
## Load modules
module load aster/1.22
# Run aster
echo "Running Aster for ${OUT_FILE_PREFIX}..."
wastral -i ${CLEAN_CAT_TREE_FILE} -o ${OUTPUT_DIR}/aster_output.tre -t ${THREAD}
if [ $? -ne 0 ]; then
    echo "Error running Aster. Please check the input files and parameters."
    exit 1
fi
echo "Aster analysis completed for ${OUT_FILE_PREFIX}. Output saved to ${OUTPUT_DIR}/aster_output.tre"
