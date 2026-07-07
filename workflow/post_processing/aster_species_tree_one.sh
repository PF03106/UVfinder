#!/bin/bash
#SBATCH --job-name=aster             # Job name
#SBATCH --mail-type=ALL              # Email notifications
#SBATCH --mail-user=                 # Replace with your email
#SBATCH --nodes=1                    # Use 1 node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=8            # Number of CPU cores
#SBATCH --mem=8gb                    # Job memory
#SBATCH --time=96:00:00              # Time limit hrs:min:sec
#SBATCH --output=logs/aster_autosomal%j.out     # Standard output and error log

# ---------------------------------------------------------------------------------------------
# 1. Set variables and working directory
OUT_FILE_PREFIX=hypnales_1kb_autosomal
THREAD=8
ANNOTATION=3

## Input directories and files
IQTREE_RESULTS_DIR1=/blue/mcdaniel/seyeonkim/UVfinder/results/hypnales_1kb/07_phylogeny/autosomal 
# ▼▼ Enter the full path to the combined metadata file below ▼▼
COMBINED_METADATA_FILE=/blue/mcdaniel/seyeonkim/UVfinder/results/metadata_dir/all_species_metadata.tsv 

## Output directory
OUTPUT_DIR=results/hypnales_1kb/${OUT_FILE_PREFIX}_aster/
mkdir -p ${OUTPUT_DIR}

## File names
BLACKLIST_FILE=${OUTPUT_DIR}/exclude_genes.txt
RAW_CAT_TREE_FILE=${OUTPUT_DIR}/${OUT_FILE_PREFIX}_raw_concatenated.treefile
CLEAN_CAT_TREE_FILE=${OUTPUT_DIR}/${OUT_FILE_PREFIX}_clean_concatenated.treefile

# -------------------------------------------------------------------------
# 2. Exclude gene trees that are duplicated, and then concatenate gene trees
echo "Step 2: Identifying duplicated genes to exclude..."

# 2-1. Create a global blacklist (using the 8th column, Gene_ID)
# [Core Logic] 
# Even if a gene (e.g., G4444) is classified as a safe 'Local_Duplication' in one species,
# if it is flagged as 'Inter_Chromosomal' (or other invalid types) in ANY other single species,
# it will be caught by the awk condition below. 
# The 'sort | uniq' command then ensures this gene is permanently excluded across all taxa.
awk -F'\t' '
$2 == "Inter_Chromosomal_Duplication" || 
$2 == "Sex_Linked_Duplication" || 
$2 == "Potential_Inter_Chrom_Dups" {print $8}
' ${COMBINED_METADATA_FILE} | sort | uniq > ${BLACKLIST_FILE}

EXCLUDED_COUNT=$(wc -l < ${BLACKLIST_FILE})
echo "Found ${EXCLUDED_COUNT} unique genes to completely exclude across all taxa. List saved to ${BLACKLIST_FILE}"

# 2-2. Concatenate only the safe gene trees that are NOT in the blacklist
echo "Concatenating clean gene trees..."
> ${RAW_CAT_TREE_FILE} # Initialize as an empty file

for tree_file in ${IQTREE_RESULTS_DIR1}/*.treefile; do
    # Extract the gene ID from the file name (e.g., G4444.treefile -> G4444)
    gene_id=$(basename ${tree_file} | grep -oE 'G[0-9]+')
    
    # If the gene_id is valid and does NOT exactly match (-qx) any line in the blacklist, append it
    if [ -n "$gene_id" ] && ! grep -qx "${gene_id}" ${BLACKLIST_FILE}; then
        cat "${tree_file}" >> ${RAW_CAT_TREE_FILE}
        echo "" >> ${RAW_CAT_TREE_FILE} # Add an empty line to prevent formatting issues between trees
    fi
done

echo "Cleanly concatenated gene trees saved to ${RAW_CAT_TREE_FILE}"
# -------------------------------------------------------------------------
# 3. Clean the concatenated tree file by removing redundant names
echo "Step 3: Cleaning sequence names in tree file..."
sed -E 's/_G[0-9]+_R[0-9]+_[AUVN]//g' ${RAW_CAT_TREE_FILE} > ${CLEAN_CAT_TREE_FILE}

# 4. Run WASTRAL(Aster)
echo "Step 4: Running Aster for ${OUT_FILE_PREFIX}..."
module load aster/1.22
wastral -i ${CLEAN_CAT_TREE_FILE} -o ${OUTPUT_DIR}/aster_output.tre -t ${THREAD} -u ${ANNOTATION}
if [ $? -ne 0 ]; then
    echo "❌ Error running Aster. Please check the input files and parameters."
    exit 1
fi

# 5. Move the output files to the output directory
if [ -f "freqQuad.csv" ]; then
    mv freqQuad.csv ${OUTPUT_DIR}/${OUT_FILE_PREFIX}_freqQuad.csv
    echo "freqQuad.csv file moved to ${OUTPUT_DIR}"
fi

echo "✔ Aster analysis completed for ${OUT_FILE_PREFIX}. Output saved to ${OUTPUT_DIR}/aster_output.tre"