#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH --array=1-334

#SBATCH --ntasks=1 # 1 core per bcf file

#SBATCH -t 1:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/07_PreFiltering

#SBATCH -e ../Job_logs/07_PreFiltering_logs/countCalls.err

#SBATCH -J countCalls

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools
module load bcftools/1.14


# Should already have a file called biallelic_bcf_paths in 07_PreFiltering directory
# File should contain the path to each bcf file (with only biallelic calls), with each path on a new line

# Pull one bcf file per task ID
# Task ID is number between 1-334. sed -n "30p" biallelic_bcf_paths would pull the path to the bcf file on line 30 of the file and store it in ${bcf}
bcf=$(sed -n "${SLURM_ARRAY_TASK_ID}p" biallelic_bcf_paths)

# Counts the number of calls (-H ignores the headers)
bcftools view -H ${bcf}| wc -l >> countCallsList # will create (initially) and append the count for each .bcf to this file

