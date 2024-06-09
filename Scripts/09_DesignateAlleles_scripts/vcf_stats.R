#!/usr/bin/env Rscript

# R script that will plot the number of variants present in each chromosome compared to chromosome length
# Example Usage:
# Rscript vcf_stats.R -t variantsPerChrom_sansMaxMiss.list -s sequence_report.tsv

# Install & load necessary packages ####

# List of packages needed
packages <- c("optparse", # v 1.7.5, allow use of command line flags
              "tidyverse") # v 2.0.0, data frame manipulation & plotting

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], repos = "https://cloud.r-project.org/") # error message if we don't set cran mirror
}

# Load the packages
invisible(lapply(packages, library, character.only = TRUE))

# Run from command line! ####

flag_list <- list(
  make_option(c("-t", "--tsv_input"), type = "character", default = NULL, 
              help = "Input TSV file that shows chromosome name (col 1) and number of variants (col 2)", metavar = "character"),
  
  make_option(c("-s", "--seqReport_input"), type = "character", default = NULL, 
              help = "Input sequence_report.tsv file, you can find it from NCBI!", metavar = "character"),
  
  make_option(c("--width"), type = "integer", default = 1600, 
              help = "Width of output images in pixels [default= %default]", metavar = "integer"),
  
  make_option(c("--height"), type = "integer", default = 1200, 
              help = "Height of output images in pixels [default= %default]", metavar = "integer"),
  
  make_option(c("--res"), type = "integer", default = 300, 
              help = "Resolution of output images in DPI [default= %default]", metavar = "integer")
)

# Set the flag list
opt_parser <- OptionParser(option_list = flag_list)
opt <- parse_args(opt_parser) # can refrence user provided arguments using opt$res for example!
  
# Read in input files ####
  
  # File with number of variants per chromosome
  variantsPerVCF <- read.csv(opt$tsv_input, header = FALSE, sep = "\t") %>%
    rename("chromosome" = "V1", 
           "number_of_variants" = "V2") 
  
  # Genome report from NCBI, https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_036417665.1/
  refGenome_report <- read.csv(opt$seqReport_input, header = TRUE, sep = "\t")

  # Genome report to include only chrom name, genbank accession, chrom length
  chromLengthDesingations_df = refGenome_report %>% 
    select(Chromosome.name, GenBank.seq.accession,Seq.length)
  
# Number of variants V chromosome length ####
  
  # Create a dataframe with chromosome name, length, and variant count all together
  nameLengthVarCount_df = full_join(variantsPerVCF, chromLengthDesingations_df, by = join_by(x$chromosome == y$GenBank.seq.accession)) %>%
    relocate(Chromosome.name, .before = chromosome) %>% # put chromosome number as 1st position
    relocate(Seq.length, .before = number_of_variants) %>% # and move sequence length to be before the site designation
    filter(chromosome != "KM078784.1",
           chromosome != "CM071464.1") # remove the mitochondrial row
    
  # Create a separate dataframe with just scaffolds
  justScaffolds_df = nameLengthVarCount_df %>%
    filter(Chromosome.name == "Un") %>% # this the name assoicated with scaffolds
    summarise(across(Seq.length,sum)) %>%
    mutate(Chromosome.name = "Un", chromosome = "scaffolds", number_of_variants = NA) %>%
    relocate(Seq.length, .before = number_of_variants)

  # concatenate all counts to create just one row repsenting scaffolds
  scaffRow = nameLengthVarCount_df %>%
    filter(chromosome == "scaffolds")
  
  # Update the dataframe to include just autosomes + Z chrom + 1 row representing scaffolds
  merged_scaffRow = scaffRow %>%
    inner_join(justScaffolds_df, by = "chromosome") %>%
    select(-Chromosome.name.x, -Seq.length.x, -number_of_variants.y) %>%
    relocate(Chromosome.name.y, .before = chromosome) %>%
    relocate(number_of_variants.x, .after = Seq.length.y) %>%
    rename("Chromosome.name" = "Chromosome.name.y",
           "Seq.length" = "Seq.length.y",
           "number_of_variants" = "number_of_variants.x"
    )
  
  nameLengthVarCount_mergedScaff_df = nameLengthVarCount_df %>%
    filter(chromosome != "scaffolds", # remove the existing scaffolds row 
           Chromosome.name != "Un") %>% # remove the separated scaffolds
    rbind(., merged_scaffRow)
  
# Plot ####
  
  png(filename = "numberOfVariants_v_chromLength.png", width = opt$width, height = opt$height, res = opt$res)
  
  ggplot(nameLengthVarCount_mergedScaff_df, aes(x = number_of_variants, y = Seq.length)) +
    geom_point(size = 3) +
    labs(title = "Number of variants in VCF vs. chrosome length",
         x = "Number of variants",
         y = "Length of chromsome") +
    scale_y_continuous(labels = scales::comma_format()) +
    scale_x_continuous(labels = scales::comma_format())
    theme_minimal()
    
  dev.off()