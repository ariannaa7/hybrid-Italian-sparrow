#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH -p core

#SBATCH -n 1

#SBATCH -t 23:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/09_DesignateAlleles/Inputs/

#SBATCH -e ../../Job_logs/09_DesignateAlleles_logs/subsetVCFbySample.err

#SBATCH -J subsetVCFsbySample

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools
module load bcftools/1.14 

for dir in ../../08_Filtering/*MaxMissing; do
    
    dirName=$(basename "$dir") # pulls just name of dir, e.g. InclMaxMissing


    for vcf_gz in ${dir}/*.vcf.gz; do

        file=$(basename "$vcf_gz" | cut -d'_' -f2) # pulls just the chromosome name or scaffolds
        
        bcftools view ${vcf_gz} --samples-file samples_to_include -Ov -o "VCFs_${dirName}/samplesubset_${file}_${dirName}.vcf"

    done
done
