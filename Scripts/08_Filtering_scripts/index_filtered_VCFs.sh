#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH -p core

#SBATCH -n 1

#SBATCH -t 15:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/08_Filtering/

#SBATCH -e ../Job_logs/08_Filtering_logs/indexing.err

#SBATCH -J indexing_filtered_VCFs

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


module load bioinfo-tools
module load bcftools/1.14


# Run samtools index on each bam file
for dir in /proj/sparrowhybridization/Pitaliae/working/Arianna/08_Filtering/*Missing*/; do # ${dir} contains entire path to directories InclMaxMissing and SansMaxMissing

    for vcf_gz in ${dir}/*.vcf.gz; do # for each file in the directory that ends in .vcf.gz

        # index each bam file and save the output in the corresponding directory with existing file name + .bai at the end (e.g. "K006_remapped.bam.bai")
        bcftools index --csi "${vcf_gz}" 

    done

done