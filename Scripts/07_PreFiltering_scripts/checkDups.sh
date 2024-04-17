#!/bin/bash

#SBATCH -A naiss2023-5-262

#SBATCH --array=1-334

#SBATCH --ntasks=1 # 1 core per bcf file

#SBATCH -t 2:00:00

#SBATCH -D /proj/sparrowhybridization/Pitaliae/working/Arianna/07_PreFiltering

#SBATCH -e ../Job_logs/07_PreFiltering_logs/checkDups.err

#SBATCH -J checkdups

#SBATCH --mail-user=ar4666al-s@student.lu.se

#SBATCH --mail-type=ALL


module load bioinfo-tools
module load bcftools/1.14

# Pull one bcf file per task ID
# Task ID is number between 1-334. sed -n "30p" biallelic_bcf_paths would pull the path to the bcf file on line 30 of the file and store it in ${bcf}
bcf=$(sed -n "${SLURM_ARRAY_TASK_ID}p" biallelic_bcf_paths)
file=$(basename "$bcf" .bcf) # will pull 126sample_JAZGKA010000334.1_biallelicSites

# check the number of calls in the bcf
beforeCheck=$(bcftools view -H "${bcf}" | cut -f 2 | wc -l) 

# check number of calls after accounting for duplciate positions
checkDup=$(bcftools view -H "${bcf}" | cut -f 2 | sort | uniq | wc -l) 

# Check if number of calls changes before and after accounting for duplicates
if [[ "$beforeCheck" != "$checkDup" ]]; # if the number of calls before and after checking for duplicates is not equal
then
    echo ${file} >> bcf_contains_dups_list # then there are duplicates present in the bcf, save the name of the bcf to this list
else
    echo "--" >> bcf_contains_dups_list  # then there aren't duplicates present in the bcf
fi