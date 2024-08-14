
# Do hybrid Italian sparrows purge derived variants disproportionately during genome stabilization?
## Overview
This repository contains the shell scripts, a Python script, and multiple R scripts that were used for exploration and to address the above question

## Data
House sparrow reference genome (BioProject PRJNA1064421), annotation file, and sequence report were accessed from NCBI
Aligned sparrow BAM files were acquired from the Runemark lab

## Workflow & scripts
One can use the workflow.mkd file to work step-by-step through the analysis, with all requirements and versions listed there. The workflow includes intermediate bash commands and refers to scripts that can be found in the associated directories nested within "Scripts". It should be noted that 11_GOanalysis_scripts was not used in the analysis due to lack of time, this is also noted in the workflow.mkd file. 

### Versions
#### Intermediate bash code found in workflow.mkd and in shell scripts requires loading the following packages:
SAMtools (v1.12), SAMtools (v1.14), BWA (v0.7.17), BCFtools (v1.14), vcftools (v0.1.16), BEDTools (v2.31.1), plink (v1.90b4.9), BEDOPS (v2.4.39), htslib (v1.14), vcflib (v1.0.1)

#### Python
Python (v3.12.1)

#### R (v4.3.3) scripts require the loading of the following packages:
optparse (v1.7.5), tidyverse (v2.0.0), eulerr (v7.0.2), stats (v4.3.3), gplots (v3.1.3.1), graphics (v4.3.3), remotes (v2.5.0), vcfR (v1.15.0), GenotypePlot (v0.2.1), ggpubr (v0.6.0), eulerr (v7.0.2), BiocManager (v1.30.23), karyoploteR (v1.28.0), dplyr (v1.1.4), tidyr (v1.3.1), magrittr (v2.0.3), rtracklayer (v1.62.0)

## Questions, comments, or concerns?
Please reach out at ar4666al-s@student.lu.se