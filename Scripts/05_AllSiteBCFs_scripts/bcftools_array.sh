#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH --array=1-334

#SBATCH --ntasks=1 # 1 core per task/region/chromosome

#SBATCH -t 4-00:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/05_AllSiteBCFs

#SBATCH -e ../Job_logs/05_AllSiteBCFs_logs/126samples_bcf.err

#SBATCH -J 126samples_bcf

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools
module load bcftools/1.14
module load htslib/1.14
module load samtools/1.14

# Save the path to the reference genome within a variable
refGenomePath=/proj/sparrowhybridization/Pitaliae/working/Arianna/Data/ReferenceGenome/GCA_036417665.1_bPasDom1.hap1_genomic.fna

# Should already have a file called mappedReads_paths and chromo_list in 05_AllSiteBCFs directory
# mappedReads_paths = contains the path to each remapped bam file on a different line
# chromo_list = contains every chromosome/region we want to generate a bcf for on a different line

# Pull one chromosome/region per task ID
# Task ID is number between 1-334. sed -n "30p" chromo_list would pull the region on line 30 of the file and store it in ${chromo}
chromo=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chromo_list)

# Run the command with the specified region
bcftools mpileup -Ou --fasta-ref ${refGenomePath} --bam-list mappedReads_paths --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR --regions ${chromo} | \
bcftools call --format-fields GQ,GP --multiallelic-caller -Ob -o 126sample_v_gc_raw_${chromo}.bcf
