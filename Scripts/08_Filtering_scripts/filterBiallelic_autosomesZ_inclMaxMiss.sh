#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH --array=1-333

#SBATCH --ntasks=20 # 20 cores per bcf file, for compression

#SBATCH -t 6:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/08_Filtering/InclMaxMissing

#SBATCH -e ../../Job_logs/08_Filtering_logs/InclMaxMissing/filtering_autosomesZ_%a.err

#SBATCH -o ../../Job_logs/08_Filtering_logs/InclMaxMissing/filtering_autosomesZ_%a.out

#SBATCH -J filtering_autosomesZ_inclMaxMiss

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools
module load vcftools/0.1.16
module load bcftools/1.14


# Pull one bcf file per task ID
# Task ID is number between 1-333. sed -n "30p" biallelic_bcf_paths_autosomesZ would pull the path to the bcf file on line 30 of the file and store it in ${bcf}
bcf=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../biallelic_bcf_paths_autosomesZ)
file=$(basename "$bcf" .bcf) # will pull 126sample_JAZGKA010000334.1_biallelicSites

# Filter
# Pipe to bcftools, following regular workflow produces a gzipped vcf that is unreadable by bcftools/vcftools
vcftools --bcf ../${bcf} --minQ 30 --min-meanDP 5 --max-meanDP 12 --minDP 3  --maxDP 18 --max-missing 0.9 \
--remove-filtered-all --recode-bcf --stdout \
| bcftools view --threads 20 -Oz -o "${file}_filtered.vcf.gz"