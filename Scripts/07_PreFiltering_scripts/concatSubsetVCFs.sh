#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH -p node

#SBATCH -t 1:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/07_PreFiltering

#SBATCH -e ../Job_logs/07_PreFiltering_logs/concat_subsets.err

#SBATCH -J concat_subsets

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


module load bioinfo-tools
module load bcftools/1.14 

# concatenate all the biallelic subset files into one compressed vcf
bcftools concat --threads 20 ./*_biallelicSites_subset.vcf -Oz -o sparrow_biallelics_subset.vcf.gz

