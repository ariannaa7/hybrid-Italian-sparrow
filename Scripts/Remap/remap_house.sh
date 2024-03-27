#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH -p node

#SBATCH -t 13:30:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/02_Remap

#SBATCH -e remapping_house.err

#SBATCH -J remap_existingBAM_house

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


module load bioinfo-tools
module load samtools/1.12 # latest version at the time of Heng Li's post
module load bwa/0.7.17 # latest version at the time of Heng Li's post


refGenomePath=/proj/sparrowhybridization/Pitaliae/working/Arianna/Data/ReferenceGenome/GCA_036417665.1_bPasDom1.hap1_genomic.fna

for bam in ../Data/Bamfiles/house/*.bam ; do # for each file in the directory that ends in .bam

    prefix=$(basename "$bam" .bam) # Pull the basename, e.g. "8M71932"

    # Code pulled from Heng Li's blogpost, see link for command-by-command description
    # https://lh3.github.io/2021/07/06/remapping-an-aligned-bam
    # 20 cores per node, 1 thread per core! Specify 20 threads for each task
    samtools collate -Oun128 "${bam}" -@ 20 | samtools fastq -OT RG,BC -@ 20 - \
        | bwa mem -pt8 -CH <(samtools view -H "${bam}"|grep ^@RG) ${refGenomePath} -t 20 - \
        | samtools sort -@ 20 -m4g -o "house_remap/${prefix}_remapped.bam" -

done