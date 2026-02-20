#!/bin/bash
#SBATCH --job-name=UVfinder           # Job name
#SBATCH --mail-type=ALL               # Email notifications
#SBATCH --mail-user=seyeonkim@ufl.edu # Replace with your email
#SBATCH --nodes=1                     # Use 1 node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=8             # Number of CPU cores
#SBATCH --mem=8gb                     # Job memory
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/RUN_all_dio_%j.out     # Standard output and error log
#SBATCH --account=mcdaniel
#SBATCH --qos=mcdaniel-b

# 1. Load Conda/Mamba
module load conda

# 2. Activate environment
# source $(conda info --base)/etc/profile.d/conda.sh
# conda activate uvfinder
conda activate /home/seyeonkim/miniconda3/envs/uvfinder

# 3. Run Snakemake
# snakemake --unlock
snakemake --cores 8 -p --latency-wait 120 --keep-going --rerun-incomplete 