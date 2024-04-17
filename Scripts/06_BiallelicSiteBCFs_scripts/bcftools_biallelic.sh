#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH --array=1-334

#SBATCH --ntasks=20 # 20 cores per bcf

#SBATCH -t 2:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/06_BiallelicSiteBCFs

#SBATCH -e ../Job_logs/06_BiallelicSiteBCFs_logs/biallelic_bcfs.err

#SBATCH -J biallelic_bcfs

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools
module load bcftools/1.14


# Should already have a file called bcf_paths in the 06_BiallelicSiteBCFs directory
# It should contain the path to each of the all site bcfs, with each path on a new lines

# Pull one bcf file per task ID
# Task ID is number between 1-334. sed -n "30p" bcf_paths would pull the path to the bcf file on line 30 of the file and store it in ${bcf}
bcf=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bcf_paths)
file=$(basename "$bcf" .bcf | sed s/_v_gc_raw//g) # will pull 126sample_v_gc_raw_JAZGKA010000334.1 then just 126sample_JAZGKA010000334.1

# pulls biallelic site, meaning 1 allele REF and 1 allele as ALT.
# -m2 = minimum 2 alleles, -M2 = max 2 alleles, -v snps = only include SNPs (extra measure to make sure indels & sites with no ALT are excluded)
bcftools view -m2 -M2 -v snps ${bcf} --threads 20 -Ob -o "${file}_biallelicSites.bcf"
