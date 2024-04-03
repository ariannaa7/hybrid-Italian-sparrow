#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH -p node

#SBATCH -t 17:30:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/02_Remap

#SBATCH -e remapping_corsica.err

#SBATCH -J remap_existingBAM_corsica

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


module load bioinfo-tools
module load samtools/1.12 # latest version at the time of Heng Li's post
module load bwa/0.7.17 # latest version at the time of Heng Li's post

# Create a subdirectory to store the output
mkdir /proj/sparrowhybridization/Pitaliae/working/Arianna/02_Remap/corsica_remap

# Save the path to the reference genome within a variable
refGenomePath=/proj/sparrowhybridization/Pitaliae/working/Arianna/Data/ReferenceGenome/GCA_036417665.1_bPasDom1.hap1_genomic.fna

for bam in ../Data/Bamfiles/corsica/*.bam ; do # for each file in the directory that ends in .bam

    prefix=$(basename "$bam" .bam | cut -d "_" -f 1) # Pull the basename then prefix, e.g. "K032_resorted_nodup_realigned" then just "K032"

    # Code pulled from Heng Li's blogpost, see link for command-by-command description
    # https://lh3.github.io/2021/07/06/remapping-an-aligned-bam
    # 20 cores per node, 1 thread per core! Specify 20 threads for each task
    samtools collate -Oun128 "${bam}" -@ 20 | samtools fastq -OT RG,BC -@ 20 - \
        | bwa mem -pt8 -CH <(samtools view -H "${bam}"|grep ^@RG) ${refGenomePath} -t 20 - \
        | samtools sort -@ 20 -m4g -o "corsica_remap/${prefix}_remapped.bam" -

done
