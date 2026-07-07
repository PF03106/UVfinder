# UVFinder: A tool to extract bryophyte sex-linked gene copies from the GoFlag408 probe set

UVfinder is a Snakemake-based bioinformatics pipeline designed to identify and map sex-linked loci in bryophyte genomes. 

## Prerequisites
Before running UVfinder, ensure you have the following installed:
* [Conda](https://docs.conda.io/en/latest/) or Mamba
* [Snakemake](https://snakemake.readthedocs.io/)

## Installation

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/PF03106/UVfinder.git](https://github.com/PF03106/UVfinder.git)
   cd UVfinder
2. **Set up the conda env**
    It is highly recommended to use Conda for environment management. You can create the required environment using the provided uvfinder.yaml
    ```bash
    conda env create -f uvfinder.yaml
3. **Configuration**
    To run UVfinder, users have to provide chromosome-level assembled genome files. Locate input genome files under 'resources/genomes' (or other preffered locations) and create samples.tsv. Below is the example of samples.tsv. (Users can put "Unknown" in the order columns)
    ```
    sample_id	order	genus	species	genome_filename
    S0001	Polytrichales	Atrichum	angustatum	GCA_050084355.1_ASM5008435v1_genomic.fna
    S0002	Hypnales	Antitrichia	curtipendula	GCA_951230915.1_cbAntCurt1.1_genomic.fna
    S0003	Andreaeales	Andreaea	nivalis	GCA_965637445.2_caAndNiva3.2_genomic.fna
    S0004	Hypnales	Aerobryopsis	subdivergens	GCA_050083985.1_ASM5008398v1_genomic.fna
    ```
    Then adjust paths(files/directories name) and parameters at config/config.yaml

4. **Running UVfinder**
    You can adjust cores, latency-wait, and other options for running UVfinder
    ```bash
        # 1. Load Conda/Mamba (if required by your institutional HPC environment)
        module load conda

        # 2. Activate the UVfinder environment
        conda activate uvfinder 
        # Note: Alternatively, use the full path to your environment (e.g., conda activate /home/username/miniconda3/envs/uvfinder)

        # 3. Unlock the working directory (if previously interrupted)
        snakemake --unlock

        # 4. Run the pipeline
        snakemake --cores 8 -p --latency-wait 120 --keep-going --rerun-incomplete
