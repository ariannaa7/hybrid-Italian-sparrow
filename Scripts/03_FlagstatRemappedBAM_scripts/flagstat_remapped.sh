#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH -p core

#SBATCH -n 16

#SBATCH -t 1:30:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/03_FlagstatRemappedBAM

#SBATCH -e ../Job_logs/03_FlagstatRemappedBAM_logs/flagstat.err

#SBATCH -J flagstat_remapped_BAM

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


module load bioinfo-tools
module load samtools/1.14

# Run samtools on each bam file
for dir in ../02_Remap/*/; do # ${dir} contains entire path

    dirName=$(basename "$dir") # just the name of the specific directory is stored (e.g corsica_remap)

    outputDir="${dirName}_flagstats" # store a directory name (e.g. corsica_remap_flagstats)

    mkdir "$outputDir" # Make a directory of this name in the current directory, 03_FlagstatRemappedBAM

    for bam in ${dir}/*.bam; do # for each file in the directory that ends in .bam

        file=$(basename "$bam" .bam); # pull the basename, e.g. "K006_remapped"

        samtools flagstat "${bam}" -@ 16 > "${outputDir}/${file}.flagstat" # run w/16 threads on each .bam file, save (e.g. "K006_remapped.flagstat") to corresponding directory

    done
done