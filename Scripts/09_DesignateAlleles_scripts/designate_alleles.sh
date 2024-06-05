#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH --array=1-41

#SBATCH --ntasks=1 # 1 core per vcf file

#SBATCH -t 2-00:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/09_DesignateAlleles/Outputs/Designations_VCFsansMaxMissing

#SBATCH -e ../../../Job_logs/09_DesignateAlleles_logs/designateAlleles.err

#SBATCH -J designateAlleles

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


# Store paths to necessary directories in environmental variables
inputsDir=/proj/sparrowhybridization/Pitaliae/working/Arianna/09_DesignateAlleles/Inputs
pythonScriptDir=/proj/sparrowhybridization/Pitaliae/working/Arianna/Scripts/09_DesignateAlleles_scripts

# Pull one vcf file per task ID
# Task ID is number between 1-41. sed -n "30p" subsetted_sansMaxMissing_VCF_paths would pull the path to the vcf file on line 30 of the file and store it in ${vcf}
vcf=$(sed -n "${SLURM_ARRAY_TASK_ID}p" subsetted_sansMaxMissing_VCF_paths)
chrom=$(basename "$vcf" | cut -d'_' -f2) # will pull just the chromosome name JAZGKA010000334



python3 ${pythonScriptDir}/sparrows_Alamshahi.py -v ${vcf} -b ${inputsDir}/popmap.tsv -c 0.80,0.85,0.90,0.95,0.99 -n ${inputsDir}/chromNames.tsv -o "${chrom}_sansMaxMiss"
