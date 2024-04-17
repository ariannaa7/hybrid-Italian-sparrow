#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH -p core

#SBATCH -n 1

#SBATCH -t 5:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/07_PreFiltering

#SBATCH -e ../Job_logs/07_PreFiltering_logs/subsetRedoZ.err

#SBATCH -J subsetRedoZ

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools
module load bcftools/1.14 
module load vcflib/1.0.1 

bcftools view ../06_BiallelicSiteBCFs/126sample_CM071456.1_biallelicSites.bcf | vcfrandomsample -r 0.025 > "Subsets/126sample_CM071456.1_biallelicSites_subset_REDO.vcf"









