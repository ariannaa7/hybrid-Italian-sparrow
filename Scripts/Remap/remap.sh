#!/bin/bash

#SBATCH -A naiss2023-5-262

#SATCH -p node

#SBATCH -t 6-20:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/02_Remap

#SBATCH -e remapping.err

#SBATCH -J remap_existingBAM

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


module load bioinfo-tools
module load samtools/1.12 # latest version at the time of Heng Li's post
module load bwa/0.7.17 # latest version at the time of Heng Li's post



refGenomePath=/proj/sparrowhybridization/Pitaliae/working/Arianna/Data/ReferenceGenome/GCA_036417665.1_bPasDom1.hap1_genomic.fna

for dir in ../Data/Bamfiles/*/; do # ${dir} contains entire path

    dirName=$(basename "$dir") # just the name of the specific directory is stored (e.g corsica)

    outputDir="${dirName}_remap" # store a directory name (e.g. corsica_remap)

    mkdir "$outputDir" # make a directory of this name in the current directory, 02_Remap

    for bam in ${dir}/*.bam; do # for each file in the directory that ends in .bam

        file=$(basename "$bam" .bam); # pull the basename, e.g. "K032_resorted_nodup_realigned"

        if [[ "$file" == Rimini* || "$file" == Lesina* ]]; # Rimini and Lesina populations follow a different format
        then
            prefix=${file} # pull the basename, e.g. "Rimini_112"
        else
            prefix=$(basename "$bam" .bam | cut -d "_" -f 1) # for all other populations, pull the basename "K032_resorted_nodup_realigned" then just "K032"
        fi

        # collate -> 
        samtools collate -Oun128 "${bam}" -@ 20 | samtools fastq -OT RG,BC -@ 20 - \
        | bwa mem -pt8 -CH <(samtools view -H "${bam}"|grep ^@RG) ${refGenomePath} -t 20 - \
        | samtools sort -@ 20 -m4g -o "${outputDir}/${prefix}_remapped.bam" -

    done
done