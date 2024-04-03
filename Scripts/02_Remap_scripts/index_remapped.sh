#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH -p core

#SBATCH -n 16

#SBATCH -t 1:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna

#SBATCH -e ../Job_logs/02_Remap_logs/indexing.err

#SBATCH -J indexing_remapped_BAM

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


module load bioinfo-tools
module load samtools/1.12 # latest version at the time of Heng Li's post


# Run samtools index on each bam file
for dir in /proj/sparrowhybridization/Pitaliae/working/Arianna/02_Remap/*/; do # ${dir} contains entire path

    for bam in ${dir}/*.bam; do # for each file in the directory that ends in .bam

        # index each bam file and save the output in the corresponding directory with existing file name + .bai at the end (e.g. "K006_remapped.bam.bai")
        samtools index "${bam}" "${bam}.bai" -@ 15 # 15 threads + main thread = 16 threads (and we are using 16 cores)

    done
done