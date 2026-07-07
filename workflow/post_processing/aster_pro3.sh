#!/bin/bash
#SBATCH --job-name=astral_pro3  # Job name
#SBATCH --mail-type=ALL                  # Email notifications
#SBATCH --mail-user=                     # Replace with your email
#SBATCH --nodes=1                        # Use 1 node
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=8                # Number of CPU cores
#SBATCH --mem=8gb                       
#SBATCH --time=96:00:00                  # Time limit hrs:min:sec
#SBATCH --output=logs/astral_pro3_%j.out   # Standard output and error log

# 1. Set variables and working directory (ADJUST PATHS IF NEEDED!)
OUT_FILE_PREFIX=sex_linked_20
ANNOTATION=3
ROOT=S0003
THREAD=$SLURM_CPUS_PER_TASK

## Input directories (iqtree results, (*.treefile))
IQTREE_RESULTS_DIR1=/blue/mcdaniel/seyeonkim/UVfinder/results_moss_new/bryo_synt/07_phylogeny/autosomal

## Output directory
OUTPUT_DIR=/blue/mcdaniel/seyeonkim/UVfinder/results_moss_new/bryo_synt/07_phylogeny/${OUT_FILE_PREFIX}_astral_pro3/
mkdir -p ${OUTPUT_DIR}

## File names
RAW_CAT_TREE_FILE=${OUTPUT_DIR}/${OUT_FILE_PREFIX}_raw_concatenated.treefile
CLEAN_CAT_TREE_FILE=${OUTPUT_DIR}/${OUT_FILE_PREFIX}_clean_concatenated.treefile

# -------------------------------------------------------------------------
# 2. Concatenate ALL gene trees (Paralogs included!)
echo "Step 1: Concatenating all gene trees..."
awk '{print}' ${IQTREE_RESULTS_DIR1}/20_collapsed*.treefile > ${RAW_CAT_TREE_FILE}
if [ $? -ne 0 ]; then
    echo "❌ Error concatenating gene trees."
    exit 1
fi
echo "All gene trees concatenated to ${RAW_CAT_TREE_FILE}"

# 3. Clean names (Crucial step for ASTRAL-pro3 multi-labelled tree format)
echo "Step 2: Cleaning sequence names for ASTRAL-pro3 auto-mapping..."
sed -E 's/_G[0-9]+_R[0-9]+_[AUVN]//g' ${RAW_CAT_TREE_FILE} > ${CLEAN_CAT_TREE_FILE}

# 4. Run ASTRAL-pro3
echo "Step 3: Running ASTRAL-Pro3..."
module load aster/1.22  

astral-pro3 -i ${CLEAN_CAT_TREE_FILE} -o ${OUTPUT_DIR}/astral_pro3_output.tre -t ${THREAD} --root ${ROOT} -u ${ANNOTATION}
if [ $? -ne 0 ]; then
    echo "❌ Error running ASTRAL-Pro3."
    exit 1
fi

echo "✔ ASTRAL-pro3 analysis completed! Output saved to ${OUTPUT_DIR}/astral_pro3_output.tre"


# 5. Move the output files to the output directory
if [ -f "freqQuad.csv" ]; then
    mv freqQuad.csv ${OUTPUT_DIR}/${OUT_FILE_PREFIX}_freqQuad.csv
    echo "freqQuad.csv file moved to ${OUTPUT_DIR}"
fi

echo "✔ Aster analysis completed for ${OUT_FILE_PREFIX}. Output saved to ${OUTPUT_DIR}/aster_output.tre"