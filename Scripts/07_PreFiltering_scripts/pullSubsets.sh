#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH --array=1-334

#SBATCH --ntasks=1 # 1 core per bcf file

#SBATCH -t 5:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/07_PreFiltering

#SBATCH -e ../Job_logs/07_PreFiltering_logs/pullSubsets.err

#SBATCH -J pullSubsets

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools
module load bcftools/1.14 
module load vcflib/1.0.1 

# Pull one bcf file per task ID
# Task ID is number between 1-334. sed -n "30p" biallelic_bcf_paths would pull the path to the bcf file on line 30 of the file and store it in ${bcf}
bcf=$(sed -n "${SLURM_ARRAY_TASK_ID}p" biallelic_bcf_paths)
file=$(basename "$bcf" .bcf) # will pull 126sample_JAZGKA010000334.1_biallelicSites

# Don't specify threads, bcftools only utilizes it for ouptut compression
# Use bcftools to view a bcf instead of vcf (goes faster than vcf or vcf.gz)
bcftools view ${bcf} | vcfrandomsample -r 0.00145 > "${file}_subset.vcf"
