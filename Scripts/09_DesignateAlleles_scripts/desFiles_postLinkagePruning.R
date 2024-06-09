#!/usr/bin/env Rscript

# Example usage! Same script as desFiles, just not using 0.80, 0.85, or 0.95 cutoffs since only using linkage pruned
# 0.99 and 0.90 cutoffs

# Rscript desFiles_postLinkagePruning.R --tenRandomCorsicans yes \
#--sansSAHA_0.90 prunedDesignationFiles/pruned_cutoff_0.90_sansSAHA.tsv \
#--SAHA_0.90 prunedDesignationFiles/pruned_cutoff_0.90_SAHA.tsv \
#--sansSAHA_0.99 prunedDesignationFiles/pruned_cutoff_0.99_sansSAHA.tsv \
#--SAHA_0.99 prunedDesignationFiles/pruned_cutoff_0.99_SAHA.tsv -s sequence_report.tsv

# Install & load necessary packages ####

# List of packages needed
packages <- c("optparse", # v 1.7.5, allow use of command line flags
              "tidyverse", # v 2.0.0, data frame manipulation & plotting
              "eulerr") # v 7.0.2, venn diagram


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], repos = "https://cloud.r-project.org/") # error message if we don't set cran mirror
}

# Load the packages
invisible(lapply(packages, library, character.only = TRUE))

# Run from command line! ####

flag_list <- list(
  
  make_option(c("--tenRandomCorsicans"), type = "character", default = NULL, 
              help = "write 'yes' or 'no'. The analysis should proceed with the 10 random corsicans = 'yes'. ", metavar = "character"),
  
  make_option(c("--sansSAHA_0.90"), type = "character", default = NULL, 
              help = "Input designation file for 0.90 cutoff, sansSAHA", metavar = "file_path"),
  
  make_option(c("--SAHA_0.90"), type = "character", default = NULL, 
              help = "Input designation file for 0.90 cutoff, SAHA", metavar = "file_path"),
  
  make_option(c("--sansSAHA_0.99"), type = "character", default = NULL, 
              help = "Input designation file for 0.99 cutoff, sansSAHA", metavar = "file_path"),
  
  make_option(c("--SAHA_0.99"), type = "character", default = NULL, 
              help = "Input designation file for 0.99 cutoff, SAHA", metavar = "file_path"),
  
  make_option(c("-s", "--seqReport_input"), type = "character", default = NULL, 
              help = "Input sequence_report.tsv file, you can find it from NCBI!", metavar = "file_path"),
  
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

# Import Data ####

# Check for random corsicans!

  #corsica_cutoff = cutoff_0.90_SAHA_dat %>% # doesn't matter 0.90 or 0.99, just using to pull names
  #select(starts_with("K0"))
  
  #ten_randomCorsicans_sample1 = sample(colnames(corsica_cutoff_resample)[1:21], 10)
  # "K039" "K083" "K073" "K034" "K026" "K022" "K031" "K010" "K085" "K019"
  ten_randomCorsicans_sample1 <- c("K039", "K083", "K073", "K034", "K026", "K022", "K031", "K010", "K085", "K019") # the 10 corsicans to include


  if (toupper(opt$tenRandomCorsicans) == "YES") {
    
    # Then the read_des function should inlcude just 10 corsicans (random)
    read_desFile <- function(file_path) {
      dat <- read.csv(file_path, header = TRUE, sep = "\t") # will replace the hash before Locus with X.
      dat <- dat %>%
        rename(
          "locus" = "X.Locus", # rename column to just "locus"
          "K034" = "X034" # rename column to "K034" (the K was dropped by vcftools when generating VCFs)
        )
      
      corsicans_to_exclude <- c("K006", "K011", "K015", "K021", "K029", "K032", "K035", "K071", "K089", "K090", "K096" ) # the 11 corsicans to exclude
      
      randomSubsetCorsicans <- dat %>%
        select(-all_of(corsicans_to_exclude))
      
      return(randomSubsetCorsicans)
    }
    
  } else if (toupper(opt$tenRandomCorsicans) == "NO") {
    
    # Then the read_des function should inlcude all the corsicans
    read_desFile <- function(file_path) {
      dat <- read.csv(file_path, header = TRUE, sep = "\t") # will replace the hash before Locus with X.
      dat <- dat %>%
        rename(
          "locus" = "X.Locus", # rename column to just "locus"
          "K034" = "X034" # rename column to "K034" (the K was dropped by vcftools when generating VCFs)
        )
      return(dat)
    }
    
    
  } else {
    stop("Looks like you didn't specify if you want to use the ten random Corsicans! Please specify yes or no")
  }

  # Function to pull just the autosomes from the designation file
  desFile_autosomes <- function(dat) {
    dat <- dat %>%
      filter(!str_detect(locus, "CM071456.1")) # remove Z chroms
    return(dat)
  }
  
  # Function to pull just the Z from the designation file
  desFile_Z <- function(dat) {
    dat <- dat %>%
      filter(str_detect(locus, "CM071456.1")) # keep only Z chrom
    return(dat)
  }

# Read in the designation files
  cutoff_0.90_sansSAHA_dat <- read_desFile(opt$sansSAHA_0.90)
  cutoff_0.90_sansSAHA_dat_autosomes <- desFile_autosomes(cutoff_0.90_sansSAHA_dat)
  cutoff_0.90_sansSAHA_dat_Z <- desFile_Z(cutoff_0.90_sansSAHA_dat)
  
  cutoff_0.90_SAHA_dat <- read_desFile(opt$SAHA_0.90)
  cutoff_0.90_SAHA_dat_autosomes <- desFile_autosomes(cutoff_0.90_SAHA_dat)
  cutoff_0.90_SAHA_dat_Z <- desFile_Z(cutoff_0.90_SAHA_dat)
  
  cutoff_0.99_sansSAHA_dat <- read_desFile(opt$sansSAHA_0.99)
  cutoff_0.99_sansSAHA_dat_autosomes <- desFile_autosomes(cutoff_0.99_sansSAHA_dat)
  cutoff_0.99_sansSAHA_dat_Z <- desFile_Z(cutoff_0.99_sansSAHA_dat)
  
  cutoff_0.99_SAHA_dat <- read_desFile(opt$SAHA_0.99)
  cutoff_0.99_SAHA_dat_autosomes <- desFile_autosomes(cutoff_0.99_SAHA_dat)
  cutoff_0.99_SAHA_dat_Z <- desFile_Z(cutoff_0.99_SAHA_dat)
 
# Genome report from NCBI, https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_036417665.1/
  refGenome_report <- read.csv(opt$seqReport_input, header = TRUE, sep = "\t")
  
# Genome report to include only chrom name, genbank accession, chrom length
  chromLengthDesingations_df <- refGenome_report %>% 
    select(Chromosome.name, GenBank.seq.accession,Seq.length)

# How many sites designated spanish derived/house ancestral & house derived/spanish ancestral at each cutoff? ####

# Function which creates a dataframes for each cutoff that shows total site SD and HD counts
  totalSite_SDHD_count <- function(dat) {
    count_df <- dat %>% # use the cutoff_[cutoff]_SAHA_dat df as base
      select(Site_Designation) %>% # select site designation column
      filter(Site_Designation %in% c("SD", "HD")) %>% # Only keep lines with SD/HD designation (precaution, but shouldn't have any -D at these cutoffs)
      group_by(Site_Designation) %>% # group by SD and HD
      summarize(count = n()) # then count total instances
    return(count_df)
  }
  
# Run the function!
  count_derivedSites_0.90 <- totalSite_SDHD_count(cutoff_0.90_SAHA_dat)
  count_derivedSites_0.90_autosomes <- totalSite_SDHD_count(cutoff_0.90_SAHA_dat_autosomes)
  count_derivedSites_0.90_Z <- totalSite_SDHD_count(cutoff_0.90_SAHA_dat_Z)
  
  count_derivedSites_0.99 <- totalSite_SDHD_count(cutoff_0.99_SAHA_dat)
  count_derivedSites_0.99_autosomes <- totalSite_SDHD_count(cutoff_0.99_SAHA_dat_autosomes)
  count_derivedSites_0.99_Z <- totalSite_SDHD_count(cutoff_0.99_SAHA_dat_Z)
  
# Create a list of all the HD/SD count dataframes to join
  count_dfs_list <- list(count_derivedSites_0.90, count_derivedSites_0.99)
  count_dfs_list_autosomes <- list(count_derivedSites_0.90_autosomes, count_derivedSites_0.99_autosomes)
  count_dfs_list_Z <- list(count_derivedSites_0.90_Z, count_derivedSites_0.99_Z)
  
  
# Initialize with the first data frame
  derivedCounts_df <- count_dfs_list[[1]] 
  derivedCounts_df_autosomes <- count_dfs_list_autosomes[[1]]
  derivedCounts_df_Z <- count_dfs_list_Z[[1]] 
  
# Join the other count dataframes to the main one with a left join
  for (i in 2:length(count_dfs_list)) {
    derivedCounts_df <- left_join(derivedCounts_df, count_dfs_list[[i]], by = "Site_Designation")
  }
  
  for (i in 2:length(count_dfs_list_autosomes)) {
    derivedCounts_df_autosomes <- left_join(derivedCounts_df_autosomes, count_dfs_list_autosomes[[i]], by = "Site_Designation")
  }
  
  for (i in 2:length(count_dfs_list_Z)) {
    derivedCounts_df_Z <- left_join(derivedCounts_df_Z, count_dfs_list_Z[[i]], by = "Site_Designation")
  }
  
# Rename the columns to match the cutoff
  derivedCounts_df <- derivedCounts_df %>%
    rename(
      "0.90" = "count.x",
      "0.99" ="count.y"
    )
  
  derivedCounts_df_autosomes <- derivedCounts_df_autosomes %>%
    rename(
      "0.90" = "count.x",
      "0.99" ="count.y"
    )
  
  derivedCounts_df_Z <- derivedCounts_df_Z %>%
    rename(
      "0.90" = "count.x",
      "0.99" ="count.y"
    )
  
# Convert to long format
  derivedCounts_df_long <- derivedCounts_df %>%
    pivot_longer(cols = -Site_Designation, names_to = "cutoff", values_to = "count")
  
  derivedCounts_df_autosomes_long <- derivedCounts_df_autosomes %>%
    pivot_longer(cols = -Site_Designation, names_to = "cutoff", values_to = "count")
  
  derivedCounts_df_Z_long <- derivedCounts_df_Z %>%
    pivot_longer(cols = -Site_Designation, names_to = "cutoff", values_to = "count")

# Functions to plot
  plot_siteDesCountsPerCutoff_bar <- function(derivedCounts_long, title) { # title should be string with ""
    plot <- ggplot(derivedCounts_long, aes(x = cutoff, y = count, fill = Site_Designation)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = "Allele frequency cutoffs", y = "Number of SNP sites with a designation", fill = "Site designation") + 
      scale_y_continuous(labels = scales::comma_format()) +
      scale_fill_manual(values = c("HD" = "skyblue1", "SD" = "sienna2"),
                        labels = c("HD" = "House derived/Spanish ancestral", "SD" = "Spanish derived/House ancestral")) +
      ggtitle(title) +
      theme_minimal() +
      theme(legend.title = element_text(hjust = 0.5)) # center the legend title

    return(plot)
  } 
  
  plot_siteDesCountsPerCutoff_stacked <- function(derivedCounts_long, title) { # title should be string with ""
    plot <- ggplot(derivedCounts_long, aes(x = cutoff, y = count, fill = Site_Designation)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(x = "Allele frequency cutoffs", y = "Number of SNP sites with a designation", fill = "Site designation") + 
      scale_y_continuous(labels = scales::comma_format()) +
      scale_fill_manual(values = c("HD" = "skyblue1", "SD" = "sienna2"),
                        labels = c("HD" = "House derived/Spanish ancestral", "SD" = "Spanish derived/House ancestral")) +
      ggtitle(title) +
      theme_minimal() +
      theme(legend.title = element_text(hjust = 0.5)) # center the legend title
    
    return(plot)
  } 
  
  # Bar Plot
  png(filename = "siteDesCountsPerCutoff_bar.png", width = opt$width, height = opt$height, res = opt$res)
  plot_siteDesCountsPerCutoff_bar(derivedCounts_df_long, "Site designation counts per cutoff")
  dev.off()
  
  png(filename = "siteDesCountsPerCutoff_autosomes_bar.png", width = opt$width, height = opt$height, res = opt$res)
  plot_siteDesCountsPerCutoff_bar(derivedCounts_df_autosomes_long, "Site designation counts per cutoff - autosomes")
  dev.off()
  
  png(filename = "siteDesCountsPerCutoff_Z_bar.png", width = opt$width, height = opt$height, res = opt$res)
  plot_siteDesCountsPerCutoff_bar(derivedCounts_df_Z_long, "Site designation counts per cutoff - Z")
  dev.off()
  
  # Stacked bar plot
  png(filename = "siteDesCountsPerCutoff_stacked.png", width = opt$width, height = opt$height, res = opt$res)
  plot_siteDesCountsPerCutoff_stacked(derivedCounts_df_long, "Site designation counts per cutoff")
  dev.off()
  
  png(filename = "siteDesCountsPerCutoff_autosomes_stacked.png", width = opt$width, height = opt$height, res = opt$res)
  plot_siteDesCountsPerCutoff_stacked(derivedCounts_df_autosomes_long, "Site designation counts per cutoff - autosomes")
  dev.off()
  
  png(filename = "siteDesCountsPerCutoff_Z_stacked.png", width = opt$width, height = opt$height, res = opt$res)
  plot_siteDesCountsPerCutoff_stacked(derivedCounts_df_Z_long, "Site designation counts per cutoff - Z")
  dev.off()
  
  # Takeaway: number of sites per decreases as cutoff increases.
  # More Spanish derived/House ancestral sites than House derived/Spanish ancestral at each site.

# How many Italian alleles designated spanish derived, house derived, or ancestral (SA or HA) at each cutoff? ####
# Create a function that will create a dataframe with total SD, HD, and A allele counts for each Italian sample
  count_Italian_alleles <- function(dat) {
    allele_count <- dat %>%
      select(7:ncol(.)) %>%
      pivot_longer(everything(), names_to = "pop", values_to = "designation") %>% # change to long format 
      group_by(pop) %>% # group by the population
      summarize(HD_count = sum(str_count(designation, "HD")),
                SD_count = sum(str_count(designation, "SD")),
                A_count = sum(str_count(designation, "SA")) + sum(str_count(designation, "HA")), # total ancestral count
                HA_count = sum(str_count(designation, "HA")), # ancestral counts split by House parentage
                SA_count = sum(str_count(designation, "SA"))) # ancestral counts split by Spanish parentage
    return(allele_count)
  }
  
# Run the function for each cutoff!
  cutoff_0.99_Ital_alleleCounts <- count_Italian_alleles(cutoff_0.99_SAHA_dat)
  cutoff_0.99_Ital_alleleCounts_autosomes <- count_Italian_alleles(cutoff_0.99_SAHA_dat_autosomes)
  cutoff_0.99_Ital_alleleCounts_Z <- count_Italian_alleles(cutoff_0.99_SAHA_dat_Z)
  
  cutoff_0.90_Ital_alleleCounts <- count_Italian_alleles(cutoff_0.90_SAHA_dat)
  cutoff_0.90_Ital_alleleCounts_autosomes <- count_Italian_alleles(cutoff_0.90_SAHA_dat_autosomes)
  cutoff_0.90_Ital_alleleCounts_Z <- count_Italian_alleles(cutoff_0.90_SAHA_dat_Z)
  
# Create a dataframe that will store total HD, SD, A counts for each cutoff
  countTotals_It_perCutoff <- data.frame() # empty dataframe
  countTotals_It_perCutoff_autosomes <- data.frame() # empty dataframe
  countTotals_It_perCutoff_Z <- data.frame() # empty dataframe
  
# List of dataframes to iterate over
  Ital_alleleCounts_dfs <- list(cutoff_0.99_Ital_alleleCounts, cutoff_0.90_Ital_alleleCounts)
  Ital_alleleCounts_dfs_autosomes <- list(cutoff_0.99_Ital_alleleCounts_autosomes, cutoff_0.90_Ital_alleleCounts_autosomes)
  Ital_alleleCounts_dfs_Z <- list(cutoff_0.99_Ital_alleleCounts_Z, cutoff_0.90_Ital_alleleCounts_Z)
  
# Loop through the list of dataframes
  for (df in Ital_alleleCounts_dfs) {
    totals <- df %>% # use the dfs we created as a base to calc totals
      summarize(HD_total = sum(HD_count), # calc the sum of each allele designation type
                SD_total = sum(SD_count),
                A_total = sum(A_count),
                HA_total = sum(HA_count),
                SA_total = sum(SA_count))
    
    # Add totals to countTotals_perCutoff
    countTotals_It_perCutoff <- bind_rows(countTotals_It_perCutoff, totals)
  }
  
  for (df in Ital_alleleCounts_dfs_autosomes) {
    totals <- df %>% # use the dfs we created as a base to calc totals
      summarize(HD_total = sum(HD_count), # calc the sum of each allele designation type
                SD_total = sum(SD_count),
                A_total = sum(A_count),
                HA_total = sum(HA_count),
                SA_total = sum(SA_count))
    
    # Add totals to countTotals_perCutoff
    countTotals_It_perCutoff_autosomes <- bind_rows(countTotals_It_perCutoff_autosomes, totals)
  }
  
  for (df in Ital_alleleCounts_dfs_Z) {
    totals <- df %>% # use the dfs we created as a base to calc totals
      summarize(HD_total = sum(HD_count), # calc the sum of each allele designation type
                SD_total = sum(SD_count),
                A_total = sum(A_count),
                HA_total = sum(HA_count),
                SA_total = sum(SA_count))
    
    # Add totals to countTotals_perCutoff
    countTotals_It_perCutoff_Z <- bind_rows(countTotals_It_perCutoff_Z, totals)
  }
  
# Append a column with cutoff values
  countTotals_It_perCutoff <- cbind(cutoff = c(0.99, 0.90), countTotals_It_perCutoff)
  countTotals_It_perCutoff_autosomes <- cbind(cutoff = c(0.99, 0.90), countTotals_It_perCutoff_autosomes)
  countTotals_It_perCutoff_Z <- cbind(cutoff = c(0.99, 0.90), countTotals_It_perCutoff_Z)
  
  
# One version of the dataframe with just ancestral, another specifying spanish ancestral and house ancestral
  countTotals_It_perCutoff_SAHA <- countTotals_It_perCutoff %>%
    select(-A_total)
  
  countTotals_It_perCutoff_sansSAHA <- countTotals_It_perCutoff %>%
    select(-HA_total, -SA_total)
  
  
  countTotals_It_perCutoff_SAHA_autosomes <- countTotals_It_perCutoff_autosomes %>%
    select(-A_total)
  
  countTotals_It_perCutoff_sansSAHA_autosomes <- countTotals_It_perCutoff_autosomes %>%
    select(-HA_total, -SA_total)
  
  
  countTotals_It_perCutoff_SAHA_Z <- countTotals_It_perCutoff_Z %>%
    select(-A_total)
  
  countTotals_It_perCutoff_sansSAHA_Z <- countTotals_It_perCutoff_Z %>%
    select(-HA_total, -SA_total)
  
# Long format!
  countTotals_It_perCutoff_SAHA_long <- pivot_longer(countTotals_It_perCutoff_SAHA, cols = -c(cutoff), names_to = "Total_type", values_to = "Total_count")
  countTotals_It_perCutoff_sansSAHA_long <- pivot_longer(countTotals_It_perCutoff_sansSAHA, cols = -c(cutoff), names_to = "Total_type", values_to = "Total_count")
  
  countTotals_It_perCutoff_SAHA_autosomes_long <- pivot_longer(countTotals_It_perCutoff_SAHA_autosomes, cols = -c(cutoff), names_to = "Total_type", values_to = "Total_count")
  countTotals_It_perCutoff_sansSAHA_autosomes_long <- pivot_longer(countTotals_It_perCutoff_sansSAHA_autosomes, cols = -c(cutoff), names_to = "Total_type", values_to = "Total_count")
  
  countTotals_It_perCutoff_SAHA_Z_long <- pivot_longer(countTotals_It_perCutoff_SAHA_Z, cols = -c(cutoff), names_to = "Total_type", values_to = "Total_count")
  countTotals_It_perCutoff_sansSAHA_Z_long <- pivot_longer(countTotals_It_perCutoff_sansSAHA_Z, cols = -c(cutoff), names_to = "Total_type", values_to = "Total_count")

  # Create plot functions
  plot_italianAlleleDesPerCutoff_sansSAHA <- function(countTotals_It_perCutoff_sansSAHA_long_df, title) { # title should be string with ""
    plot <-  ggplot(countTotals_It_perCutoff_sansSAHA_long_df, aes(x = as.factor(cutoff), y = Total_count, fill = Total_type)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = "Allele frequency cutoffs", y = "Total number of Italian alleles", fill = "Italian allele designation") +
      scale_y_continuous(labels = scales::comma_format()) +
      scale_fill_manual(values = c("HD_total" = "royalblue", "SD_total" = "orangered3", "A_total" = "orchid"),
                        labels = c("HD_total" = "House derived", "SD_total" = "Spanish derived", "A_total" = "Ancestral")) +
      scale_x_discrete(labels = c("0.9" = "0.90")) + # rename these cutoffs to have two decimal places
      ggtitle(title) +
      theme_minimal()
    
    return(plot)
  } 
  
  plot_italianAlleleDesPerCutoff_SAHA <- function(countTotals_It_perCutoff_SAHA_long_df, title) { # title should be string with ""
    plot <-  ggplot(countTotals_It_perCutoff_SAHA_long_df, aes(x = as.factor(cutoff), y = Total_count, fill = Total_type)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = "Allele frequency cutoffs", y = "Total number of Italian alleles", fill = "Italian allele designation") +
      scale_y_continuous(labels = scales::comma_format()) +
      scale_fill_manual(values = c("HD_total" = "royalblue", "SD_total" = "orangered3", "HA_total" = "#5c2be2", "SA_total" = "magenta2"),
                        labels = c("HD_total" = "House derived", "SD_total" = "Spanish derived", "HA_total" = "House ancestral", "SA_total" = "Spanish ancestral")) +
      scale_x_discrete(labels = c("0.9" = "0.90")) + # rename these cutoffs to have two decimal places
      ggtitle(title) +
      theme_minimal()
    
    return(plot)
  } 
  
  # Plot sansSAHA
  png(filename = "italianAlleleDesPerCutoff_sansSAHA.png", width = opt$width, height = opt$height, res = opt$res)
  plot_italianAlleleDesPerCutoff_sansSAHA(countTotals_It_perCutoff_sansSAHA_long, "Total allele counts in Italians per cutoff" )
  dev.off()
  
  png(filename = "italianAlleleDesPerCutoff_sansSAHA_autosomes.png", width = opt$width, height = opt$height, res = opt$res)
  plot_italianAlleleDesPerCutoff_sansSAHA(countTotals_It_perCutoff_sansSAHA_autosomes_long, "Total allele counts in Italians per cutoff - autosomes" )
  dev.off()
  
  png(filename = "italianAlleleDesPerCutoff_sansSAHA_Z.png", width = opt$width, height = opt$height, res = opt$res)
  plot_italianAlleleDesPerCutoff_sansSAHA(countTotals_It_perCutoff_sansSAHA_Z_long, "Total allele counts in Italians per cutoff - Z" )
  dev.off()
  
  
# Plot SAHA
  png(filename = "italianAlleleDesPerCutoff_SAHA.png", width = opt$width, height = opt$height, res = opt$res)
  plot_italianAlleleDesPerCutoff_SAHA(countTotals_It_perCutoff_SAHA_long, "Total allele counts in Italians per cutoff")
  dev.off()
  
  png(filename = "italianAlleleDesPerCutoff_SAHA_autosomes.png", width = opt$width, height = opt$height, res = opt$res)
  plot_italianAlleleDesPerCutoff_SAHA(countTotals_It_perCutoff_SAHA_autosomes_long, "Total allele counts in Italians per cutoff - autosomes")
  dev.off()
  
  png(filename = "italianAlleleDesPerCutoff_SAHA_Z.png", width = opt$width, height = opt$height, res = opt$res)
  plot_italianAlleleDesPerCutoff_SAHA(countTotals_It_perCutoff_SAHA_Z_long, "Total allele counts in Italians per cutoff - Z")
  dev.off()


  
  # Takeaways: use 0.99 as hard cutoff and 0.90 as soft cutoff. 
  # More HD than SD at higher cutoff, flips at 0.90. Not a reason to choose 0.90, not that interesting
  # Always have more ancestral alleles, specifically more house ancestral
  
# How many Italian alleles are designated having spanish or house parentage at each cutoff? ####
  
  # Function that creates a dataframe with total S and h allele counts for each Italian sample
  justSH_countItalian <- function(allele_counts_df, cutoffValue) { # enter the cutoffValue as a string in ""
    SH_count <- allele_counts_df %>%
      mutate(S_count = rowSums(select(., contains("S"))),
             H_count = rowSums(select(., contains("H"))),
             cutoff = cutoffValue) %>%
      select(-HA_count, -HD_count, -SA_count, -SD_count, -A_count) %>%
      relocate(cutoff, .before = pop)
    
    return(SH_count)
  }
  
  # Run the function for each cutoff!
  cutoff_0.90_Ital_alleleCounts_justSH = justSH_countItalian(cutoff_0.90_Ital_alleleCounts, "0.90")
  cutoff_0.90_Ital_alleleCounts_justSH_autosomes = justSH_countItalian(cutoff_0.90_Ital_alleleCounts_autosomes, "0.90")
  cutoff_0.90_Ital_alleleCounts_justSH_Z = justSH_countItalian(cutoff_0.90_Ital_alleleCounts_Z, "0.90")
  
  cutoff_0.99_Ital_alleleCounts_justSH = justSH_countItalian(cutoff_0.99_Ital_alleleCounts, "0.99")
  cutoff_0.99_Ital_alleleCounts_justSH_autosomes = justSH_countItalian(cutoff_0.99_Ital_alleleCounts_autosomes, "0.99")
  cutoff_0.99_Ital_alleleCounts_justSH_Z = justSH_countItalian(cutoff_0.99_Ital_alleleCounts_Z, "0.99")
  
  # Create a function that converts each of the above dataframes to long format
  justSH_countItalian_long <- function(justSH_countItalian_df) {
    long_format <- justSH_countItalian_df %>%
      pivot_longer(cols = starts_with("S_count") | starts_with("H_count"), names_to = "parent_count", values_to = "count")
    return(long_format)
  }
  
  # Run the function for each cutoff!
  cutoff_0.90_Ital_alleleCounts_justSH_long <- justSH_countItalian_long(cutoff_0.90_Ital_alleleCounts_justSH)
  cutoff_0.90_Ital_alleleCounts_justSH_autosomes_long <- justSH_countItalian_long(cutoff_0.90_Ital_alleleCounts_justSH_autosomes)
  cutoff_0.90_Ital_alleleCounts_justSH_Z_long <- justSH_countItalian_long(cutoff_0.90_Ital_alleleCounts_justSH_Z)
  
  cutoff_0.99_Ital_alleleCounts_justSH_long <- justSH_countItalian_long(cutoff_0.99_Ital_alleleCounts_justSH)
  cutoff_0.99_Ital_alleleCounts_justSH_autosomes_long <- justSH_countItalian_long(cutoff_0.99_Ital_alleleCounts_justSH_autosomes)
  cutoff_0.99_Ital_alleleCounts_justSH_Z_long <- justSH_countItalian_long(cutoff_0.99_Ital_alleleCounts_justSH_Z)
  
  
  # Combine into one df!
  Ital_alleleCounts_justSH_long <- bind_rows(cutoff_0.90_Ital_alleleCounts_justSH_long,
                                            cutoff_0.99_Ital_alleleCounts_justSH_long)
  
  Ital_alleleCounts_justSH_autosomes_long <- bind_rows(cutoff_0.90_Ital_alleleCounts_justSH_autosomes_long,
                                            cutoff_0.99_Ital_alleleCounts_justSH_autosomes_long)
  
  Ital_alleleCounts_justSH_Z_long <- bind_rows(cutoff_0.90_Ital_alleleCounts_justSH_Z_long,
                                            cutoff_0.99_Ital_alleleCounts_justSH_Z_long)
  
  # Plot function
  plot_italAlleleParentCountPerCutoff <- function(Ital_alleleCounts_justSH_long_df, title) { # title should be string with ""
    plot <- ggplot(Ital_alleleCounts_justSH_long_df, aes(x = as.factor(cutoff), y = count, fill = parent_count)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = "Allele frequency cutoffs", y = "Total number of Italian alleles", fill = "Italian allele parentage") +
      scale_y_continuous(labels = scales::comma_format()) +
      scale_fill_manual(values = c("H_count" = "#FF7B65", "S_count" = "#9BD184"),
                        labels = c("H_count" = "House", "S_count" = "Spanish")) +
      scale_x_discrete(labels = c("0.9" = "0.90")) + # rename these cutoffs to have two decimal places
      ggtitle(title) +
      theme_minimal()
    
    return(plot)
  } 
  
  # Plot
  png(filename = "italAlleleParentCountPerCutoff.png", width = opt$width, height = opt$height, res = opt$res)
  plot_italAlleleParentCountPerCutoff(Ital_alleleCounts_justSH_long, "Italian allele parentage counts per cutoff")
  dev.off()
  
  png(filename = "italAlleleParentCountPerCutoff_autosomes.png", width = opt$width, height = opt$height, res = opt$res)
  plot_italAlleleParentCountPerCutoff(Ital_alleleCounts_justSH_autosomes_long, 'Italian allele parentage counts per cutoff - autosomes')
  dev.off()
  
  png(filename = "italAlleleParentCountPerCutoff_Z.png", width = opt$width, height = opt$height, res = opt$res)
  plot_italAlleleParentCountPerCutoff(Ital_alleleCounts_justSH_Z_long, 'Italian allele parentage counts per cutoff - Z')
  dev.off()
  
  #Takeaways: overall, more house alleles present
  
# For each italian population how many alleles designated spanish derived, house derived, or ancestral @ 0.90 & 0.99 ####

# Function to create a dataframe that groups samples by population prefix and shows total alleles designated HD, SD, HA, SA, or A
  ital_prefix_allele_counts <- function(allele_counts_df) {
    ital_prefix_allele_counts_df <- allele_counts_df %>%
      mutate(pop_prefix = substr(pop, 1, 2)) %>% # pull the first two letters from pop, that is your population prefix
      relocate(pop_prefix, .before = pop) %>% # move the prefix column before the pop column
      select(-pop) %>% # remove pop column
      group_by(pop_prefix) %>% # group by pop prefix (so, by population)
      summarize(
        HD_total = sum(HD_count), # then pull total counts for each designation
        SD_total = sum(SD_count),
        A_total = sum(A_count),
        HA_total = sum(HA_count),
        SA_total = sum(SA_count)
      )
    return(ital_prefix_allele_counts_df)
  }
  
# Run the function for strict (0.99) and loose (0.90) cutoff
  cutoff_0.99_ItalPrefix_alleleCounts <- ital_prefix_allele_counts(cutoff_0.99_Ital_alleleCounts)
  cutoff_0.99_ItalPrefix_alleleCounts_autosomes <- ital_prefix_allele_counts(cutoff_0.99_Ital_alleleCounts_autosomes)
  cutoff_0.99_ItalPrefix_alleleCounts_Z <- ital_prefix_allele_counts(cutoff_0.99_Ital_alleleCounts_Z)
  
  cutoff_0.90_ItalPrefix_alleleCounts <- ital_prefix_allele_counts(cutoff_0.90_Ital_alleleCounts)
  cutoff_0.90_ItalPrefix_alleleCounts_autosomes <- ital_prefix_allele_counts(cutoff_0.90_Ital_alleleCounts_autosomes)
  cutoff_0.90_ItalPrefix_alleleCounts_Z <- ital_prefix_allele_counts(cutoff_0.90_Ital_alleleCounts_Z)
  
# Create two functions, one to pull just ancestral in df and another one to specify SA and HA in df
  ital_prefix_allele_counts_SAHA <- function(ItalPrefix_alleleCounts_df) {
    ItalPrefix_alleleCounts_SAHA <- ItalPrefix_alleleCounts_df %>%
      select(-A_total)
    return(ItalPrefix_alleleCounts_SAHA)
  }
  
  ital_prefix_allele_counts_sansSAHA <- function(ItalPrefix_alleleCounts_df) {
    ItalPrefix_alleleCounts_sansSAHA <- ItalPrefix_alleleCounts_df %>%
      select(-HA_total, -SA_total)
    return(ItalPrefix_alleleCounts_sansSAHA)
  }
  
# Run the functions!
  cutoff_0.99_ItalPrefix_alleleCounts_SAHA <- ital_prefix_allele_counts_SAHA(cutoff_0.99_ItalPrefix_alleleCounts)
  cutoff_0.99_ItalPrefix_alleleCounts_SAHA_autosomes <- ital_prefix_allele_counts_SAHA(cutoff_0.99_ItalPrefix_alleleCounts_autosomes)
  cutoff_0.99_ItalPrefix_alleleCounts_SAHA_Z <- ital_prefix_allele_counts_SAHA(cutoff_0.99_ItalPrefix_alleleCounts_Z)
  
  cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA <- ital_prefix_allele_counts_sansSAHA(cutoff_0.99_ItalPrefix_alleleCounts)
  cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA_autosomes <- ital_prefix_allele_counts_sansSAHA(cutoff_0.99_ItalPrefix_alleleCounts_autosomes)
  cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA_Z <- ital_prefix_allele_counts_sansSAHA(cutoff_0.99_ItalPrefix_alleleCounts_Z)
  
  cutoff_0.90_ItalPrefix_alleleCounts_SAHA <- ital_prefix_allele_counts_SAHA(cutoff_0.90_ItalPrefix_alleleCounts)
  cutoff_0.90_ItalPrefix_alleleCounts_SAHA_autosomes <- ital_prefix_allele_counts_SAHA(cutoff_0.90_ItalPrefix_alleleCounts_autosomes) 
  cutoff_0.90_ItalPrefix_alleleCounts_SAHA_Z <- ital_prefix_allele_counts_SAHA(cutoff_0.90_ItalPrefix_alleleCounts_Z) 
  
  cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA <- ital_prefix_allele_counts_sansSAHA(cutoff_0.90_ItalPrefix_alleleCounts)
  cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA_autosomes <- ital_prefix_allele_counts_sansSAHA(cutoff_0.90_ItalPrefix_alleleCounts_autosomes)
  cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA_Z <- ital_prefix_allele_counts_sansSAHA(cutoff_0.90_ItalPrefix_alleleCounts_Z)

# Create a function to turn each of the four dataframes into long format!
  ItalPrefix_alleleCounts_long <- function(ItalPrefix_alleleCounts_SAHA_orSansSAHA_df) {
    long_format <- ItalPrefix_alleleCounts_SAHA_orSansSAHA_df %>%
      pivot_longer(cols = -c(pop_prefix), names_to = "Total_type", values_to = "Total_count")
    return(long_format)
  }
  
# Run the function
  cutoff_0.99_ItalPrefix_alleleCounts_SAHA_long <- ItalPrefix_alleleCounts_long(cutoff_0.99_ItalPrefix_alleleCounts_SAHA)
  cutoff_0.99_ItalPrefix_alleleCounts_SAHA_autosomes_long <- ItalPrefix_alleleCounts_long(cutoff_0.99_ItalPrefix_alleleCounts_SAHA_autosomes)
  cutoff_0.99_ItalPrefix_alleleCounts_SAHA_Z_long <- ItalPrefix_alleleCounts_long(cutoff_0.99_ItalPrefix_alleleCounts_SAHA_Z)
  
  cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA_long <- ItalPrefix_alleleCounts_long(cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA)
  cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA_autosomes_long <- ItalPrefix_alleleCounts_long(cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA_autosomes)
  cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA_Z_long <- ItalPrefix_alleleCounts_long(cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA_Z)
  
  cutoff_0.90_ItalPrefix_alleleCounts_SAHA_long <- ItalPrefix_alleleCounts_long(cutoff_0.90_ItalPrefix_alleleCounts_SAHA)
  cutoff_0.90_ItalPrefix_alleleCounts_SAHA_autosomes_long <- ItalPrefix_alleleCounts_long(cutoff_0.90_ItalPrefix_alleleCounts_SAHA_autosomes)
  cutoff_0.90_ItalPrefix_alleleCounts_SAHA_Z_long <- ItalPrefix_alleleCounts_long(cutoff_0.90_ItalPrefix_alleleCounts_SAHA_Z)
  
  cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA_long <- ItalPrefix_alleleCounts_long(cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA)
  cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA_autosomes_long <- ItalPrefix_alleleCounts_long(cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA_autosomes)
  cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA_Z_long <- ItalPrefix_alleleCounts_long(cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA_Z)
  
  
# Create two functions to make a stacked bar plot (one for sansSAHA, one for SAHA) 
  plot_ItalPrefix_alleleCounts_sansSAHA <- function(ItalPrefix_alleleCounts_sansSAHA_long, title) { # title should be string with ""
    plot <- ggplot(ItalPrefix_alleleCounts_sansSAHA_long, aes(x = as.factor(pop_prefix), y = Total_count, fill = Total_type)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(x = "Italian population", y = "Total number of alleles", fill = "Italian allele designation") +
      scale_y_continuous(labels = scales::comma_format()) +
      scale_fill_manual(values = c("HD_total" = "royalblue", "SD_total" = "orangered3", "A_total" = "orchid"),
                        labels = c("HD_total" = "House derived", "SD_total" = "Spanish derived", "A_total" = "Ancestral")) +
      scale_x_discrete(labels = c("C0" = "Crete", "K0" = "Corsica", "M0" = "Malta", "S0" = "Sicily" )) +
      ggtitle(title) +
      theme_minimal()
    return(plot)
  } 
  
  plot_ItalPrefix_alleleCounts_SAHA <- function(ItalPrefix_alleleCounts_SAHA_long, title) { # title should be string with ""
    plot <- ggplot(ItalPrefix_alleleCounts_SAHA_long, aes(x = as.factor(pop_prefix), y = Total_count, fill = Total_type)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(x = "Italian population", y = "Total number of alleles", fill = "Italian allele designation") +
      scale_y_continuous(labels = scales::comma_format()) +
      scale_fill_manual(values = c("HD_total" = "royalblue", "SD_total" = "orangered3", "HA_total" = "#5c2be2", "SA_total" = "magenta2"),
                        labels = c("HD_total" = "House derived", "SD_total" = "Spanish derived", "HA_total" = "House ancestral", "SA_total" = "Spanish ancestral")) +
      scale_x_discrete(labels = c("C0" = "Crete", "K0" = "Corsica", "M0" = "Malta", "S0" = "Sicily" )) +
      ggtitle(title) +
      theme_minimal()
    return(plot)
  } 
  
  # Plot!
  
  png(filename = "italAlleleDesPerPop_0.90_sansSAHA.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_sansSAHA(cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA_long, title = "Total allele counts per Italian population at 0.90 cutoff")
  dev.off()
  
  png(filename = "italAlleleDesPerPop_0.90_sansSAHA_autosomes.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_sansSAHA(cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA_autosomes_long, title = "Total allele counts per Italian population at 0.90 cutoff - autosomes")
  dev.off()
  
  png(filename = "italAlleleDesPerPop_0.90_sansSAHA_Z.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_sansSAHA(cutoff_0.90_ItalPrefix_alleleCounts_sansSAHA_Z_long, title = "Total allele counts per Italian population at 0.90 cutoff - Z")
  dev.off()
  
  
  
  png(filename = "italAlleleDesPerPop_0.99_sansSAHA.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_sansSAHA(cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA_long, title = "Total allele counts per Italian population at 0.99 cutoff")
  dev.off()
  
  png(filename = "italAlleleDesPerPop_0.99_sansSAHA_autosomes.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_sansSAHA(cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA_autosomes_long, title = "Total allele counts per Italian population at 0.99 cutoff - autosomes")
  dev.off()
  
  png(filename = "italAlleleDesPerPop_0.99_sansSAHA_Z.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_sansSAHA(cutoff_0.99_ItalPrefix_alleleCounts_sansSAHA_Z_long, title = "Total allele counts per Italian population at 0.99 cutoff - Z")
  dev.off()
  
  
  
  png(filename = "italAlleleDesPerPop_0.90_SAHA.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_SAHA(cutoff_0.90_ItalPrefix_alleleCounts_SAHA_long, title = "Total allele counts per Italian population at 0.90 cutoff")
  dev.off()
  
  png(filename = "italAlleleDesPerPop_0.90_SAHA_autosomes.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_SAHA(cutoff_0.90_ItalPrefix_alleleCounts_SAHA_autosomes_long, title = "Total allele counts per Italian population at 0.90 cutoff - autosomes")
  dev.off()
  
  png(filename = "italAlleleDesPerPop_0.90_SAHA_Z.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_SAHA(cutoff_0.90_ItalPrefix_alleleCounts_SAHA_Z_long, title = "Total allele counts per Italian population at 0.90 cutoff - Z")
  dev.off()
  
  
  
  png(filename = "italAlleleDesPerPop_0.99_SAHA.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_SAHA(cutoff_0.99_ItalPrefix_alleleCounts_SAHA_long, title = "Total allele counts per Italian population at 0.99 cutoff")
  dev.off()
  
  png(filename = "italAlleleDesPerPop_0.99_SAHA_autosomes.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_SAHA(cutoff_0.99_ItalPrefix_alleleCounts_SAHA_autosomes_long, title = "Total allele counts per Italian population at 0.99 cutoff - autosomes")
  dev.off()
  
  png(filename = "italAlleleDesPerPop_0.99_SAHA_Z.png", width = opt$width, height = opt$height, res = opt$res)
  plot_ItalPrefix_alleleCounts_SAHA(cutoff_0.99_ItalPrefix_alleleCounts_SAHA_Z_long, title = "Total allele counts per Italian population at 0.99 cutoff - Z")
  dev.off()
  
  # Takeways: Majority of Italian alleles are Ancestral and of those ancestral, most are house
  # Crete and corsica we see more house alleles, both derived and ancestral > matches PCA/admixture
  # Malta & sicily we see more spanish alleles, both derived and ancestral > matches PCA/admixture
  
# For each italian population, how many are heterozygote/homozygote at 0.90/0/.99 cutoff ####

# Create a function that will make dataframe with all the italian genotypes in the designation file grouped by population
  italGenotypes <- function(dat) {
    ital_genotypes_df <- dat %>%
      select(7:ncol(.)) %>%
      pivot_longer(cols = everything(), names_to = "sample", values_to = "genotype") %>%
      mutate(pop_prefix = substr(sample, 1, 2)) %>%
      relocate(pop_prefix, .before = sample) %>%
      select(-sample) %>%
      group_by(pop_prefix) %>%
      mutate(rowId = row_number()) %>%  # need to add row ID or else pivot wider won't work because there are duplicate genotypes and population names
      pivot_wider(names_from = pop_prefix, values_from = genotype) %>%
      select(-rowId) # will produce NAs because more K0 samples than all others
    
    return(ital_genotypes_df)
  }
  
# Run the function for 0.90/0.99 cutoffs SAHA/sansSAHA
  cutoff_0.90_Ital_genotypes_SAHA <- italGenotypes(cutoff_0.90_SAHA_dat)
  cutoff_0.90_Ital_genotypes_SAHA_autosomes <- italGenotypes(cutoff_0.90_SAHA_dat_autosomes)
  cutoff_0.90_Ital_genotypes_SAHA_Z <- italGenotypes(cutoff_0.90_SAHA_dat_Z)
  
  cutoff_0.90_Ital_genotypes_sansSAHA <- italGenotypes(cutoff_0.90_sansSAHA_dat)
  cutoff_0.90_Ital_genotypes_sansSAHA_autosomes <- italGenotypes(cutoff_0.90_sansSAHA_dat_autosomes)
  cutoff_0.90_Ital_genotypes_sansSAHA_Z <- italGenotypes(cutoff_0.90_sansSAHA_dat_Z)
  
  cutoff_0.99_Ital_genotypes_SAHA <- italGenotypes(cutoff_0.99_SAHA_dat)
  cutoff_0.99_Ital_genotypes_SAHA_autosomes <- italGenotypes(cutoff_0.99_SAHA_dat_autosomes)
  cutoff_0.99_Ital_genotypes_SAHA_Z <- italGenotypes(cutoff_0.99_SAHA_dat_Z)
  
  cutoff_0.99_Ital_genotypes_sansSAHA <- italGenotypes(cutoff_0.99_sansSAHA_dat)
  cutoff_0.99_Ital_genotypes_sansSAHA_autosomes <- italGenotypes(cutoff_0.99_sansSAHA_dat_autosomes)
  cutoff_0.99_Ital_genotypes_sansSAHA_Z <- italGenotypes(cutoff_0.99_sansSAHA_dat_Z)
  
# Function to create a dataframe with all the counts of each possible genotype grouped by population
  italGenotypeCounts <- function(italGenotypes_df) {
    italGenotypesCounts_df <- italGenotypes_df %>%
      pivot_longer(cols = everything(), names_to = "population", values_to = "genotype") %>% # transform to long format
      group_by(population, genotype) %>% # group by population and genotype
      summarize(count = n()) %>% # then count occurrence of each genotype for each population
      pivot_wider(names_from = genotype, values_from = count, values_fill = list(count = 0)) %>% # change back to wide format to easily remove NA column
      select(-"NA") # remove the NA column, only exists because more K0 samples than all others
    
    return(italGenotypesCounts_df)
  }
  
# Run the function for each cutoff & SAHA/sansSAHA combo
  genotype_counts_0.90_SAHA <- italGenotypeCounts(cutoff_0.90_Ital_genotypes_SAHA)
  genotype_counts_0.90_SAHA_autosomes <- italGenotypeCounts(cutoff_0.90_Ital_genotypes_SAHA_autosomes)
  genotype_counts_0.90_SAHA_Z <- italGenotypeCounts(cutoff_0.90_Ital_genotypes_SAHA_Z)
  
  genotype_counts_0.90_sansSAHA <- italGenotypeCounts(cutoff_0.90_Ital_genotypes_sansSAHA)
  genotype_counts_0.90_sansSAHA_autosomes <- italGenotypeCounts(cutoff_0.90_Ital_genotypes_sansSAHA_autosomes)
  genotype_counts_0.90_sansSAHA_Z <- italGenotypeCounts(cutoff_0.90_Ital_genotypes_sansSAHA_Z)
  
  genotype_counts_0.99_SAHA <- italGenotypeCounts(cutoff_0.99_Ital_genotypes_SAHA)
  genotype_counts_0.99_SAHA_autosomes <- italGenotypeCounts(cutoff_0.99_Ital_genotypes_SAHA_autosomes)
  genotype_counts_0.99_SAHA_Z <- italGenotypeCounts(cutoff_0.99_Ital_genotypes_SAHA_Z)
  
  genotype_counts_0.99_sansSAHA <- italGenotypeCounts(cutoff_0.99_Ital_genotypes_sansSAHA)
  genotype_counts_0.99_sansSAHA_autosomes <- italGenotypeCounts(cutoff_0.99_Ital_genotypes_sansSAHA_autosomes)
  genotype_counts_0.99_sansSAHA_Z <- italGenotypeCounts(cutoff_0.99_Ital_genotypes_sansSAHA_Z)
    
# Functions to create a dataframe with all the counts of each possible type of zygosity grouped by population
  zygosityCounts_sansSAHA <- function(genotype_counts_sansSAHA_df) {
    zygosity_counts_sansSAHA_df <- genotype_counts_sansSAHA_df %>%
      mutate(heterozygous_houseDerived_spanishAncestral = `A/HD` + `HD/A`)%>%
      mutate(heterozygous_spanishDerived_houseAncestral = `A/SD` + `SD/A`) %>% # if genotype in VCF 0/1 and in site des alt was derived or if ref was derived
      rename(
        "missing_genotype" = "-", # rename the genotypes to long name
        "homozygous_ancestral" = "A/A",
        "homozygous_houseDerived" = "HD/HD",
        "homozygous_spanishDerived" = "SD/SD"
      ) %>%
      select(-`A/SD`,-`SD/A`, -`A/HD`, -`HD/A`) %>% # remove the genotype columns we used for our calculations
      mutate(genotype_total = missing_genotype + homozygous_ancestral + # Add a total column at the end
               homozygous_houseDerived + homozygous_spanishDerived + 
               heterozygous_houseDerived_spanishAncestral + heterozygous_spanishDerived_houseAncestral)
    
    return(zygosity_counts_sansSAHA_df)
  }
  
  zygosityCounts_SAHA <- function(genotype_counts_SAHA_df) {
    zygosity_counts_SAHA_df <- genotype_counts_SAHA_df %>%
      mutate(heterozygous_houseDerived_spanishAncestral = `HD/SA` + `SA/HD`)%>%
      mutate(heterozygous_spanishDerived_houseAncestral = `SD/HA` + `HA/SD`) %>% # if genotype in VCF 0/1 and in site des alt was derived or if ref was derived
      rename(
        "missing_genotype" = "-", # rename the genotypes to long name
        "homozygous_houseAncestral" = "HA/HA",
        "homozygous_spanishAncestral" = "SA/SA",
        "homozygous_houseDerived" = "HD/HD",
        "homozygous_spanishDerived" = "SD/SD"
      ) %>%
      select(-`HA/SD`,-`SD/HA`, -`SA/HD`, -`HD/SA`) %>%
      mutate(genotype_total = missing_genotype + homozygous_houseAncestral + # Add a total column at the end
               homozygous_spanishAncestral + homozygous_houseDerived + 
               homozygous_spanishDerived + heterozygous_houseDerived_spanishAncestral + 
               heterozygous_spanishDerived_houseAncestral)
    
    return(zygosity_counts_SAHA_df)
  }
  
# Run the functions for each cutoff & SAHA/sansSAHA combo
  zygosity_counts_0.90_sansSAHA <- zygosityCounts_sansSAHA(genotype_counts_0.90_sansSAHA)
  zygosity_counts_0.90_sansSAHA_autosomes <- zygosityCounts_sansSAHA(genotype_counts_0.90_sansSAHA_autosomes)
  zygosity_counts_0.90_sansSAHA_Z <- zygosityCounts_sansSAHA(genotype_counts_0.90_sansSAHA_Z)
  
  zygosity_counts_0.99_sansSAHA <- zygosityCounts_sansSAHA(genotype_counts_0.99_sansSAHA)
  zygosity_counts_0.99_sansSAHA_autosomes <- zygosityCounts_sansSAHA(genotype_counts_0.99_sansSAHA_autosomes)
  zygosity_counts_0.99_sansSAHA_Z <- zygosityCounts_sansSAHA(genotype_counts_0.99_sansSAHA_Z)
  
  zygosity_counts_0.90_SAHA <- zygosityCounts_SAHA(genotype_counts_0.90_SAHA)
  zygosity_counts_0.90_SAHA_autosomes <- zygosityCounts_SAHA(genotype_counts_0.90_SAHA_autosomes)
  zygosity_counts_0.90_SAHA_Z <- zygosityCounts_SAHA(genotype_counts_0.90_SAHA_Z)
  
  zygosity_counts_0.99_SAHA <- zygosityCounts_SAHA(genotype_counts_0.99_SAHA)
  zygosity_counts_0.99_SAHA_autosomes <- zygosityCounts_SAHA(genotype_counts_0.99_SAHA_autosomes)
  zygosity_counts_0.99_SAHA_Z <- zygosityCounts_SAHA(genotype_counts_0.99_SAHA_Z)
  
# Functions to create a dataframe of frequency per zygosity type per population
  zygosity_freq_sansSAHA <- function (zygosity_counts_sansSAHA_df) {
    zygosity_freq_sansSAHA_df <- zygosity_counts_sansSAHA_df %>%
      rowwise() %>%
      mutate( # add columns that show frequency
        freq_missing_genotype = missing_genotype / genotype_total,
        freq_homozygous_ancestral = homozygous_ancestral / genotype_total,
        freq_homozygous_houseDerived = homozygous_houseDerived / genotype_total,
        freq_homozygous_spanishDerived = homozygous_spanishDerived / genotype_total,
        freq_heterozygous_houseDerived_spanishAncestral = heterozygous_houseDerived_spanishAncestral / genotype_total,
        freq_heterozygous_spanishDerived_houseAncestral = heterozygous_spanishDerived_houseAncestral / genotype_total
      ) %>%
      # remove the count columns, leaving just population and frequency columns
      select(-missing_genotype, -homozygous_ancestral, -homozygous_houseDerived, -homozygous_spanishDerived, -heterozygous_houseDerived_spanishAncestral, -heterozygous_spanishDerived_houseAncestral, -genotype_total) %>% 
      ungroup() # ungroup data to create stacked bar chart!
    
    return(zygosity_freq_sansSAHA_df)
  }
  
  zygosity_freq_SAHA <- function (zygosity_counts_SAHA_df) {
    zygosity_freq_SAHA_df <- zygosity_counts_SAHA_df %>%
      rowwise() %>%
      mutate( # add columns that show frequency
        freq_missing_genotype = missing_genotype / genotype_total,
        freq_homozygous_houseAncestral = homozygous_houseAncestral / genotype_total,
        freq_homozygous_spanishAncestral = homozygous_spanishAncestral / genotype_total,
        freq_homozygous_houseDerived = homozygous_houseDerived / genotype_total,
        freq_homozygous_spanishDerived = homozygous_spanishDerived / genotype_total,
        freq_heterozygous_houseDerived_spanishAncestral = heterozygous_houseDerived_spanishAncestral / genotype_total,
        freq_heterozygous_spanishDerived_houseAncestral = heterozygous_spanishDerived_houseAncestral / genotype_total
      ) %>%
      # remove the count columns, leaving just population and frequency columns
      select(-missing_genotype, -homozygous_houseAncestral, -homozygous_spanishAncestral, -homozygous_houseDerived, -homozygous_spanishDerived, -heterozygous_houseDerived_spanishAncestral, -heterozygous_spanishDerived_houseAncestral, -genotype_total) %>% 
      ungroup() # ungroup data to create stacked bar chart!
    
    return(zygosity_freq_SAHA_df)
  }
  
# Run the functions on each 0.90/0.99 SAHA/sansSAHA combo
  zygosity_freq_0.90_sansSAHA <- zygosity_freq_sansSAHA(zygosity_counts_0.90_sansSAHA)
  zygosity_freq_0.90_sansSAHA_autosomes <- zygosity_freq_sansSAHA(zygosity_counts_0.90_sansSAHA_autosomes)
  zygosity_freq_0.90_sansSAHA_Z <- zygosity_freq_sansSAHA(zygosity_counts_0.90_sansSAHA_Z)
  
  zygosity_freq_0.99_sansSAHA <- zygosity_freq_sansSAHA(zygosity_counts_0.99_sansSAHA)
  zygosity_freq_0.99_sansSAHA_autosomes <- zygosity_freq_sansSAHA(zygosity_counts_0.99_sansSAHA_autosomes)
  zygosity_freq_0.99_sansSAHA_Z <- zygosity_freq_sansSAHA(zygosity_counts_0.99_sansSAHA_Z)
  
  zygosity_freq_0.90_SAHA <- zygosity_freq_SAHA(zygosity_counts_0.90_SAHA)
  zygosity_freq_0.90_SAHA_autosomes <- zygosity_freq_SAHA(zygosity_counts_0.90_SAHA_autosomes)
  zygosity_freq_0.90_SAHA_Z <- zygosity_freq_SAHA(zygosity_counts_0.90_SAHA_Z)
  
  zygosity_freq_0.99_SAHA <- zygosity_freq_SAHA(zygosity_counts_0.99_SAHA)
  zygosity_freq_0.99_SAHA_autosomes <- zygosity_freq_SAHA(zygosity_counts_0.99_SAHA_autosomes)
  zygosity_freq_0.99_SAHA_Z <- zygosity_freq_SAHA(zygosity_counts_0.99_SAHA_Z)
  
  
# Create a function to convert to long format!
  zygosity_freq_long <- function(zygosity_freq_SAHA_orSansSAHA_df) {
    long_format <- zygosity_freq_SAHA_orSansSAHA_df %>%
      pivot_longer(cols = starts_with("freq"), names_to = "genotype", values_to = "frequency")
    return(long_format)
  }

# Run the function on each 0.90/0.99 SAHA/sansSAHA combo
  zygosity_freq_0.90_sansSAHA_long <- zygosity_freq_long(zygosity_freq_0.90_sansSAHA)
  zygosity_freq_0.90_sansSAHA_autosomes_long <- zygosity_freq_long(zygosity_freq_0.90_sansSAHA_autosomes)
  zygosity_freq_0.90_sansSAHA_Z_long <- zygosity_freq_long(zygosity_freq_0.90_sansSAHA_Z)
  
  zygosity_freq_0.99_sansSAHA_long <- zygosity_freq_long(zygosity_freq_0.99_sansSAHA)
  zygosity_freq_0.99_sansSAHA_autosomes_long <- zygosity_freq_long(zygosity_freq_0.99_sansSAHA_autosomes)
  zygosity_freq_0.99_sansSAHA_Z_long <- zygosity_freq_long(zygosity_freq_0.99_sansSAHA_Z)
  
  zygosity_freq_0.90_SAHA_long <- zygosity_freq_long(zygosity_freq_0.90_SAHA)
  zygosity_freq_0.90_SAHA_autosomes_long <- zygosity_freq_long(zygosity_freq_0.90_SAHA_autosomes)
  zygosity_freq_0.90_SAHA_Z_long <- zygosity_freq_long(zygosity_freq_0.90_SAHA_Z)
  
  zygosity_freq_0.99_SAHA_long <- zygosity_freq_long(zygosity_freq_0.99_SAHA)
  zygosity_freq_0.99_SAHA_autosomes_long <- zygosity_freq_long(zygosity_freq_0.99_SAHA_autosomes)
  zygosity_freq_0.99_SAHA_Z_long <- zygosity_freq_long(zygosity_freq_0.99_SAHA_Z)
  
  
  
# Create functions to plot stacked bar chart (SAHA/sansSAHA)
  plot_zygosity_freq_sansSAHA <- function(zygosity_freq_sansSAHA_long, title) { # title should be string with ""
    plot <- ggplot(zygosity_freq_sansSAHA_long, aes(x = population, y = frequency, fill = genotype)) +
      geom_bar(stat = "identity") +
      labs(x = "Italian population", y = "Genotype frequency", fill = "Genotype") +
      scale_fill_manual(values = c("freq_missing_genotype" = "lawngreen", "freq_homozygous_ancestral" = "orchid", "freq_homozygous_houseDerived" = "royalblue", "freq_homozygous_spanishDerived" = "orangered3", "freq_heterozygous_houseDerived_spanishAncestral" = "skyblue1", "freq_heterozygous_spanishDerived_houseAncestral" = "sienna2"), 
                        labels = c("freq_missing_genotype" = "Missing genotype", "freq_homozygous_ancestral" = "Homozygous: Ancestral", "freq_homozygous_houseDerived" = "Homozgyous: House Derived", "freq_homozygous_spanishDerived" = "Homozygous: Spanish Derived", "freq_heterozygous_houseDerived_spanishAncestral" = "Heterozygous: House Derived/Spanish Ancestral", "freq_heterozygous_spanishDerived_houseAncestral" = "Heterozygous: Spanish Derived/House Ancestral")) +
      scale_x_discrete(labels = c("C0" = "Crete", "K0" = "Corsica", "M0" = "Malta", "S0" = "Sicily" )) +
      ggtitle(title) +
      theme_minimal()
    return(plot)
  } 
  
  plot_zygosity_freq_SAHA <- function(zygosity_freq_SAHA_long, title) { # title should be string with ""
    plot <- ggplot(zygosity_freq_SAHA_long, aes(x = population, y = frequency, fill = genotype)) +
      geom_bar(stat = "identity") +
      labs(x = "Italian population", y = "Genotype frequency", fill = "Genotype") +
      scale_fill_manual(values = c("freq_missing_genotype" = "lawngreen", "freq_homozygous_houseAncestral" = "#5c2be2", "freq_homozygous_spanishAncestral" = "magenta2", "freq_homozygous_houseDerived" = "royalblue", "freq_homozygous_spanishDerived" = "orangered3", "freq_heterozygous_houseDerived_spanishAncestral" = "skyblue1", "freq_heterozygous_spanishDerived_houseAncestral" = "sienna2"), 
                        labels = c("freq_missing_genotype" = "Missing genotype", "freq_homozygous_houseAncestral" = "Homozygous: House Ancestral", "freq_homozygous_spanishAncestral" = "Homozygous: Spanish Ancestral", "freq_homozygous_houseDerived" = "Homozgyous: House Derived", "freq_homozygous_spanishDerived" = "Homozygous: Spanish Derived", "freq_heterozygous_houseDerived_spanishAncestral" = "Heterozygous: House Derived/Spanish Ancestral", "freq_heterozygous_spanishDerived_houseAncestral" = "Heterozygous: Spanish Derived/House Ancestral")) +
      scale_x_discrete(labels = c("C0" = "Crete", "K0" = "Corsica", "M0" = "Malta", "S0" = "Sicily" )) +
      ggtitle(title) +
      theme_minimal()
    return(plot)
  } 
  
# Plot the charts using the functions
  
  png(filename = "genotypeFreqPerPop_0.90_sansSAHA.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_sansSAHA(zygosity_freq_0.90_sansSAHA_long, "Genotype frequency per population at 0.90 cutoff")
  dev.off()
  
  png(filename = "genotypeFreqPerPop_0.90_sansSAHA_autosomes.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_sansSAHA(zygosity_freq_0.90_sansSAHA_autosomes_long, "Genotype frequency per population at 0.90 cutoff - autosomes")
  dev.off()
  
  png(filename = "genotypeFreqPerPop_0.90_sansSAHA_Z.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_sansSAHA(zygosity_freq_0.90_sansSAHA_Z_long, "Genotype frequency per population at 0.90 cutoff - Z")
  dev.off()
  
  
  
  png(filename = "genotypeFreqPerPop_0.99_sansSAHA.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_sansSAHA(zygosity_freq_0.99_sansSAHA_long, "Genotype frequency per population at 0.99 cutoff")
  dev.off()
  
  png(filename = "genotypeFreqPerPop_0.99_sansSAHA_autosomes.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_sansSAHA(zygosity_freq_0.99_sansSAHA_autosomes_long, "Genotype frequency per population at 0.99 cutoff - autosomes")
  dev.off()
  
  png(filename = "genotypeFreqPerPop_0.99_sansSAHA_Z.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_sansSAHA(zygosity_freq_0.99_sansSAHA_Z_long, "Genotype frequency per population at 0.99 cutoff - Z")
  dev.off()
  
  
  
  png(filename = "genotypeFreqPerPop_0.90_SAHA.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_SAHA(zygosity_freq_0.90_SAHA_long, "Genotype frequency per population at 0.90 cutoff")
  dev.off()
  
  png(filename = "genotypeFreqPerPop_0.90_SAHA_autosomes.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_SAHA(zygosity_freq_0.90_SAHA_autosomes_long, "Genotype frequency per population at 0.90 cutoff - autosomes")
  dev.off()
  
  png(filename = "genotypeFreqPerPop_0.90_SAHA_Z.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_SAHA(zygosity_freq_0.90_SAHA_Z_long, "Genotype frequency per population at 0.90 cutoff - Z")
  dev.off()
  
  
  
  png(filename = "genotypeFreqPerPop_0.99_SAHA.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_SAHA(zygosity_freq_0.99_SAHA_long, "Genotype frequency per population at 0.99 cutoff")
  dev.off()
  
  png(filename = "genotypeFreqPerPop_0.99_SAHA_autosomes.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_SAHA(zygosity_freq_0.99_SAHA_autosomes_long, "Genotype frequency per population at 0.99 cutoff - autosomes")
  dev.off()
  
  png(filename = "genotypeFreqPerPop_0.99_SAHA_Z.png", width = 2000, height = 1300, res = opt$res)
  plot_zygosity_freq_SAHA(zygosity_freq_0.99_SAHA_Z_long, "Genotype frequency per population at 0.99 cutoff - Z")
  dev.off()
  
  # Takeaways: similar frequency of homozygous ancestral across populations.
  # Crete and corsica had more homozygous house derived > matches PCA/admixture 
  # Crete and corsica more heterozygous house ancestral than other het
  # Malta and sicily had more homozygous spanish derived > matches PCA/admixture
  # Malta and sicily more heterozygous spanish derived than other het

# Number of types of sites v. chrom length @ 0.90 & 0.99 cutoff ####

# Keep code the way it was originally where I removed Z chroms from full designation set!
# Function to subset cutoff_0.90_sansSAHA_dat & cutoff_0.99_sansSAHA_dat to include chrom names & site designation (no allele freqs or italian designations)
  nameDes_subset <- function(dat) {
    nameDes_subset_df <- dat %>% 
      select(locus, Site_Designation) %>%
      mutate(chromosome = sapply(strsplit(as.character(locus), "_"), `[`, 1)) %>% # add column with just chromname (genbank)
      relocate(chromosome, .before = locus) %>% # move the genbank ID chrom column to first position
      select(-locus) # don't actually need locus, just needed it to pull genbank accession/ID
    return(nameDes_subset_df)
  }
  
# Run the function!
  nameDes_0.90_subset <- nameDes_subset(cutoff_0.90_sansSAHA_dat)
  nameDes_0.99_subset <- nameDes_subset(cutoff_0.99_sansSAHA_dat)
  
# Create functions to: get associated seq length, join the refgenome subsetted dataframe and nameDes_0.90_subset/nameDes_0.99_subset dataframe based on genbank accession
  nameDesLength <- function(nameDes_subset,chromLengthDesingations){
    nameDes_subset = full_join(nameDes_subset, chromLengthDesingations, by = join_by(x$chromosome == y$GenBank.seq.accession)) %>%
      relocate(Chromosome.name, .before = chromosome) %>% # put chromosome number as 1st position
      relocate(Seq.length, .before = Site_Designation) # and move sequence length to be before the site designation
    return(nameDes_subset) # there will be more rows than present in designation file bc scaffolds with designated sites (& mt) are included
  }
  
  desVlength <- function(nameDesLength_df) {
    desVlength_DF = nameDesLength_df %>%
      group_by(chromosome) %>% # group based on genbank accession chrom ID
      summarize(
        Seq.length = unique(Seq.length), # Each chrom is listed multiple times (based on locus), but length is sames so just pull once
        SD_count = sum(Site_Designation == "SD"), # count number of sites designated spanish derived
        HD_count = sum(Site_Designation == "HD") # count number of sites designated house derived
      ) %>%
      arrange(Seq.length) %>% # sort based on chrom length, ascending
      filter(!chromosome == "CM071456.1") %>% # remove sex chromosome Z, only want to look at autsomes for this (we know that Z will have a LOT of derived sites!)
      na.omit() # remove rows without site designations (short scaffolds etc.)
    return(desVlength_DF) # CHECK exclamation point placement for filter
  }

# Run the functions
  nameDesLength_0.90 <- nameDesLength(nameDes_0.90_subset,chromLengthDesingations_df) # chrom number, name, length of chrom, & site designation
  nameDesLength_0.99 <- nameDesLength(nameDes_0.99_subset,chromLengthDesingations_df)
  
  desVlength_0.90_DF <- desVlength(nameDesLength_0.90) # chrom name, chrom length, # of SD sites on chrom, # of HD sites on chrom
  desVlength_0.99_DF <- desVlength(nameDesLength_0.99)
  
# Function to convert to long format based on SD_Count/HD_Count
  desVlength_long <- function(desVlength_DF) {
    long_format <- desVlength_DF %>%
      pivot_longer(cols = c(SD_count, HD_count), names_to = "Designation", values_to = "Count")
    return(long_format)
  }
  
# Run function to convert to long format
  desVlength_0.90_DF_long <- desVlength_long(desVlength_0.90_DF)
  desVlength_0.99_DF_long <- desVlength_long(desVlength_0.99_DF)
  
  
# Create a plot function
  plot_desVlength <- function(desVlength_DF_long, title) { # title should be string with ""
    plot <- ggplot(desVlength_DF_long, aes(x = Count, y = Seq.length, color = Designation)) +
      geom_point(size = 3) +
      labs(x = "Number of sites",
           y = "Length of chromsome",
           color = "Site designation type") +
      scale_y_continuous(labels = scales::comma_format()) +
      scale_color_manual(values = c("SD_count" = "skyblue1", "HD_count" = "sienna2"),
                         labels = c("HD_count" = "House derived/Spanish Ancestral", "SD_count" = "Spanish derived/House ancestral")) +
      ggtitle(title) +
      theme_minimal()
    return(plot)
  } 

# Run the plot function for the two types of cutoffs
  
  png(filename = "sitesvChromLength_0.90_SDHD.png", width = 2000, height = 1300, res = opt$res)
  plot_desVlength(desVlength_0.90_DF_long, "Number of sites vs. chromosome length (0.90 cutoff - autosomes)")
  dev.off()
  
  png(filename = "sitesvChromLength_0.99_SDHD.png", width = 2000, height = 1300, res = opt$res)
  plot_desVlength(desVlength_0.99_DF_long, "Number of sites vs. chromosome length (0.99 cutoff - autosomes)")
  dev.off()
  
# Number of sites v. chrom length @ 0.90 & 0.99 cutoff ####
  
# Create function to find the total number of designations (regardless of parentage) per chromosome
  
  desGenVlength_DF <- function(desVlength_DF) {
    desGenVlength_DF <- desVlength_DF %>%
      mutate(total_designations = SD_count + HD_count) %>%
      select(-SD_count,-HD_count)
  }
  
# Run the function!
  desGenVlength_0.90_DF <- desGenVlength_DF(desVlength_0.90_DF)
  desGenVlength_0.99_DF <- desGenVlength_DF(desVlength_0.99_DF)
  
# Create a plot function
  plot_desGenVlength <- function(desGenVlengt_DF, title) { # title should be string with ""
    plot <- ggplot(desGenVlengt_DF, aes(x = total_designations, y = Seq.length)) +
      geom_point(size = 3) +
      labs(x = "Number of sites",
           y = "Length of chromsome") +
      scale_y_continuous(labels = scales::comma_format()) +
      ggtitle(title) +
      theme_minimal()
    return(plot)
  } 

# Plot!
  png(filename = "sitesvChromLength_0.90.png", width = 2000, height = 1300, res = opt$res)
  plot_desGenVlength(desGenVlength_0.90_DF, "Number of sites vs. chromosome length (0.90 cutoff - autosomes)")
  dev.off()

  png(filename = "sitesvChromLength_0.99.png", width = 2000, height = 1300, res = opt$res)
  plot_desGenVlength(desGenVlength_0.99_DF, "Number of sites vs. chromosome length (0.99 cutoff - autosomes)")
  dev.off()
  
# Venn diagrams! ####
  
# For each population (at 0.90 & 0.99 cutoff), create a list with locus and allele (split up genotype) for every sample
  vennList <- function(dat, popPrefix) { # popPrefix should be string in ""
    vennList_df <- dat %>%
      select(locus, starts_with(popPrefix)) %>%
      gather(key = popPrefix, value = "genotype", -locus) %>% # long format, pulls only specified population and lists every locus and genotype
      separate(genotype, into = c("allele1", "allele2"), sep = "/") %>% # separates genotype into two alleles
      gather(key = "allele_Pos", value = "allele", allele1, allele2) %>% # shows locus, sample, allele position 1 or 2 (since we split genotype), and the allele
      select(locus, allele) %>% # just keep the locus and the allele
      unique() # keep unique instance (e.g homozygous SD/SD, then would have two entries with same locus/allele pair)
    return(vennList_df)
  }
  
  corsica_vennList_cutoff_0.90 <- vennList(cutoff_0.90_SAHA_dat, popPrefix = "K0")
  corsica_vennList_cutoff_0.90_autosomes <- vennList(cutoff_0.90_SAHA_dat_autosomes, popPrefix = "K0")
  corsica_vennList_cutoff_0.90_Z <- vennList(cutoff_0.90_SAHA_dat_Z, popPrefix = "K0")
  
  crete_vennList_cutoff_0.90 <- vennList(cutoff_0.90_SAHA_dat, popPrefix = "C0")
  crete_vennList_cutoff_0.90_autosomes <- vennList(cutoff_0.90_SAHA_dat_autosomes, popPrefix = "C0")
  crete_vennList_cutoff_0.90_Z <- vennList(cutoff_0.90_SAHA_dat_Z, popPrefix = "C0")
  
  malta_vennList_cutoff_0.90 <- vennList(cutoff_0.90_SAHA_dat, popPrefix = "M0")
  malta_vennList_cutoff_0.90_autosomes <- vennList(cutoff_0.90_SAHA_dat_autosomes, popPrefix = "M0")
  malta_vennList_cutoff_0.90_Z <- vennList(cutoff_0.90_SAHA_dat_Z, popPrefix = "M0")
  
  sicily_vennList_cutoff_0.90 <- vennList(cutoff_0.90_SAHA_dat, popPrefix = "S0")
  sicily_vennList_cutoff_0.90_autosomes <- vennList(cutoff_0.90_SAHA_dat_autosomes, popPrefix = "S0")
  sicily_vennList_cutoff_0.90_Z <- vennList(cutoff_0.90_SAHA_dat_Z, popPrefix = "S0")
  
  
  corsica_vennList_cutoff_0.99 <- vennList(cutoff_0.99_SAHA_dat, popPrefix = "K0")
  corsica_vennList_cutoff_0.99_autosomes <- vennList(cutoff_0.99_SAHA_dat_autosomes, popPrefix = "K0")
  corsica_vennList_cutoff_0.99_Z <- vennList(cutoff_0.99_SAHA_dat_Z, popPrefix = "K0")
  
  crete_vennList_cutoff_0.99 <- vennList(cutoff_0.99_SAHA_dat, popPrefix = "C0")
  crete_vennList_cutoff_0.99_autosomes <- vennList(cutoff_0.99_SAHA_dat_autosomes, popPrefix = "C0")
  crete_vennList_cutoff_0.99_Z <- vennList(cutoff_0.99_SAHA_dat_Z, popPrefix = "C0")
  
  malta_vennList_cutoff_0.99 <- vennList(cutoff_0.99_SAHA_dat, popPrefix = "M0")
  malta_vennList_cutoff_0.99_autosomes <- vennList(cutoff_0.99_SAHA_dat_autosomes, popPrefix = "M0")
  malta_vennList_cutoff_0.99_Z <- vennList(cutoff_0.99_SAHA_dat_Z, popPrefix = "M0")
  
  sicily_vennList_cutoff_0.99 <- vennList(cutoff_0.99_SAHA_dat, popPrefix = "S0")
  sicily_vennList_cutoff_0.99_autosomes <- vennList(cutoff_0.99_SAHA_dat_autosomes, popPrefix = "S0")
  sicily_vennList_cutoff_0.99_Z <- vennList(cutoff_0.99_SAHA_dat_Z, popPrefix = "S0")
  

# Use eulerr package to create a fit list
  
  # Look at any occurence of SD at the locus (can be heterozygote/homozygote)
  fit_list_SD_0.90 <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.90$locus[corsica_vennList_cutoff_0.90$allele == "SD"]),
                                       Crete = unique(crete_vennList_cutoff_0.90$locus[crete_vennList_cutoff_0.90$allele == "SD"]),
                                       Malta = unique(malta_vennList_cutoff_0.90$locus[malta_vennList_cutoff_0.90$allele == "SD"]),
                                       Sicily = unique(sicily_vennList_cutoff_0.90$locus[sicily_vennList_cutoff_0.90$allele == "SD"]))
                                  )
  
  fit_list_SD_0.90_autosomes <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.90_autosomes$locus[corsica_vennList_cutoff_0.90_autosomes$allele == "SD"]),
                                        Crete = unique(crete_vennList_cutoff_0.90_autosomes$locus[crete_vennList_cutoff_0.90_autosomes$allele == "SD"]),
                                        Malta = unique(malta_vennList_cutoff_0.90_autosomes$locus[malta_vennList_cutoff_0.90_autosomes$allele == "SD"]),
                                        Sicily = unique(sicily_vennList_cutoff_0.90_autosomes$locus[sicily_vennList_cutoff_0.90_autosomes$allele == "SD"]))
  )
  
  fit_list_SD_0.90_Z <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.90_Z$locus[corsica_vennList_cutoff_0.90_Z$allele == "SD"]),
                                        Crete = unique(crete_vennList_cutoff_0.90_Z$locus[crete_vennList_cutoff_0.90_Z$allele == "SD"]),
                                        Malta = unique(malta_vennList_cutoff_0.90_Z$locus[malta_vennList_cutoff_0.90_Z$allele == "SD"]),
                                        Sicily = unique(sicily_vennList_cutoff_0.90_Z$locus[sicily_vennList_cutoff_0.90_Z$allele == "SD"]))
  )
  
  
  
  fit_list_SD_0.99 <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.99$locus[corsica_vennList_cutoff_0.99$allele == "SD"]),
                                       Crete = unique(crete_vennList_cutoff_0.99$locus[crete_vennList_cutoff_0.99$allele == "SD"]),
                                       Malta = unique(malta_vennList_cutoff_0.99$locus[malta_vennList_cutoff_0.99$allele == "SD"]),
                                       Sicily = unique(sicily_vennList_cutoff_0.99$locus[sicily_vennList_cutoff_0.99$allele == "SD"]))
                                   )
  
  fit_list_SD_0.99_autosomes <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.99_autosomes$locus[corsica_vennList_cutoff_0.99_autosomes$allele == "SD"]),
                                        Crete = unique(crete_vennList_cutoff_0.99_autosomes$locus[crete_vennList_cutoff_0.99_autosomes$allele == "SD"]),
                                        Malta = unique(malta_vennList_cutoff_0.99_autosomes$locus[malta_vennList_cutoff_0.99_autosomes$allele == "SD"]),
                                        Sicily = unique(sicily_vennList_cutoff_0.99_autosomes$locus[sicily_vennList_cutoff_0.99_autosomes$allele == "SD"]))
  )
  
  fit_list_SD_0.99_Z <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.99_Z$locus[corsica_vennList_cutoff_0.99_Z$allele == "SD"]),
                                        Crete = unique(crete_vennList_cutoff_0.99_Z$locus[crete_vennList_cutoff_0.99_Z$allele == "SD"]),
                                        Malta = unique(malta_vennList_cutoff_0.99_Z$locus[malta_vennList_cutoff_0.99_Z$allele == "SD"]),
                                        Sicily = unique(sicily_vennList_cutoff_0.99_Z$locus[sicily_vennList_cutoff_0.99_Z$allele == "SD"]))
  )
  
  # Look at any occurence of HD at the locus (can be heterozygote/homozygote)
  fit_list_HD_0.90 <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.90$locus[corsica_vennList_cutoff_0.90$allele == "HD"]),
                                       Crete = unique(crete_vennList_cutoff_0.90$locus[crete_vennList_cutoff_0.90$allele == "HD"]),
                                       Malta = unique(malta_vennList_cutoff_0.90$locus[malta_vennList_cutoff_0.90$allele == "HD"]),
                                       Sicily = unique(sicily_vennList_cutoff_0.90$locus[sicily_vennList_cutoff_0.90$allele == "HD"]))
                                   )
  
  fit_list_HD_0.90_autosomes <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.90_autosomes$locus[corsica_vennList_cutoff_0.90_autosomes$allele == "HD"]),
                                        Crete = unique(crete_vennList_cutoff_0.90_autosomes$locus[crete_vennList_cutoff_0.90_autosomes$allele == "HD"]),
                                        Malta = unique(malta_vennList_cutoff_0.90_autosomes$locus[malta_vennList_cutoff_0.90_autosomes$allele == "HD"]),
                                        Sicily = unique(sicily_vennList_cutoff_0.90_autosomes$locus[sicily_vennList_cutoff_0.90_autosomes$allele == "HD"]))
  )
  
  fit_list_HD_0.90_Z <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.90_Z$locus[corsica_vennList_cutoff_0.90_Z$allele == "HD"]),
                                        Crete = unique(crete_vennList_cutoff_0.90_Z$locus[crete_vennList_cutoff_0.90_Z$allele == "HD"]),
                                        Malta = unique(malta_vennList_cutoff_0.90_Z$locus[malta_vennList_cutoff_0.90_Z$allele == "HD"]),
                                        Sicily = unique(sicily_vennList_cutoff_0.90_Z$locus[sicily_vennList_cutoff_0.90_Z$allele == "HD"]))
  )
  
  
  
  fit_list_HD_0.99 <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.99$locus[corsica_vennList_cutoff_0.99$allele == "HD"]),
                                        Crete = unique(crete_vennList_cutoff_0.99$locus[crete_vennList_cutoff_0.99$allele == "HD"]),
                                        Malta = unique(malta_vennList_cutoff_0.99$locus[malta_vennList_cutoff_0.99$allele == "HD"]),
                                        Sicily = unique(sicily_vennList_cutoff_0.99$locus[sicily_vennList_cutoff_0.99$allele == "HD"]))
                                    )
  
  fit_list_HD_0.99_autosomes <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.99_autosomes$locus[corsica_vennList_cutoff_0.99_autosomes$allele == "HD"]),
                                        Crete = unique(crete_vennList_cutoff_0.99_autosomes$locus[crete_vennList_cutoff_0.99_autosomes$allele == "HD"]),
                                        Malta = unique(malta_vennList_cutoff_0.99_autosomes$locus[malta_vennList_cutoff_0.99_autosomes$allele == "HD"]),
                                        Sicily = unique(sicily_vennList_cutoff_0.99_autosomes$locus[sicily_vennList_cutoff_0.99_autosomes$allele == "HD"]))
  )
  
  fit_list_HD_0.99_Z <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.99_Z$locus[corsica_vennList_cutoff_0.99_Z$allele == "HD"]),
                                        Crete = unique(crete_vennList_cutoff_0.99_Z$locus[crete_vennList_cutoff_0.99_Z$allele == "HD"]),
                                        Malta = unique(malta_vennList_cutoff_0.99_Z$locus[malta_vennList_cutoff_0.99_Z$allele == "HD"]),
                                        Sicily = unique(sicily_vennList_cutoff_0.99_Z$locus[sicily_vennList_cutoff_0.99_Z$allele == "HD"]))
  )
  
# Create plot function for venn diagram!
  plot_venn <- function(fitList, title) { # title should be string with ""
    plot <- plot(fitList, # provide your fit list you made above
                 quantities = TRUE, # the plot will include the #s in each part of the diagram
                 fill = c("Corsica" = "#E9F5F9", "Crete" = "#9FBFD5", "Malta" = "#ED9668", "Sicily" = "#FAE2A0"),
                 main = title )
    return(plot)
  } 
  
# Plot
  png(filename = "vennSD_0.90.png", width = 1700, height = 2100, res = opt$res)
  plot_venn(fit_list_SD_0.90, "Spanish derived allele present in at least one copy\nat the same locus across samples in populations\n(0.90 cutoff)")
  dev.off()
  
  png(filename = "vennSD_0.90_autosomes.png", width = 1700, height = 2100, res = opt$res)
  plot_venn(fit_list_SD_0.90_autosomes, "Spanish derived allele present in at least one copy\nat the same locus across samples in populations\n(0.90 cutoff - autosomes)")
  dev.off()
  
  png(filename = "vennSD_0.90_Z.png", width = 1700, height = 2100, res = opt$res)
  plot_venn(fit_list_SD_0.90_Z, "Spanish derived allele present in at least one copy\nat the same locus across samples in populations\n(0.90 cutoff - Z)")
  dev.off()
  
  
  
  png(filename = "vennHD_0.90.png", width = 1650, height = 2100, res = opt$res)
  plot_venn(fit_list_HD_0.90, "House derived allele present in at least one copy\nat the same locus across samples in populations\n(0.90 cutoff)")
  dev.off()
  
  png(filename = "vennHD_0.90_autosomes.png", width = 1650, height = 2100, res = opt$res)
  plot_venn(fit_list_HD_0.90_autosomes, "House derived allele present in at least one copy\nat the same locus across samples in populations\n(0.90 cutoff - autosomes)")
  dev.off()
  
  png(filename = "vennHD_0.90_Z.png", width = 1650, height = 2100, res = opt$res)
  plot_venn(fit_list_HD_0.90_Z, "House derived allele present in at least one copy\nat the same locus across samples in populations\n(0.90 cutoff - Z)")
  dev.off()
  
  
  
  png(filename = "vennSD_0.99.png", width = 1700, height = 2100, res = opt$res)
  plot_venn(fit_list_SD_0.99, "Spanish derived allele present in at least one copy\nat the same locus across samples in populations\n(0.99 cutoff)")
  dev.off()
  
  png(filename = "vennSD_0.99_autosomes.png", width = 1700, height = 2100, res = opt$res)
  plot_venn(fit_list_SD_0.99_autosomes, "Spanish derived allele present in at least one copy\nat the same locus across samples in populations\n(0.99 cutoff - autosomes)")
  dev.off()
  
  png(filename = "vennSD_0.99_Z.png", width = 1700, height = 2100, res = opt$res)
  plot_venn(fit_list_SD_0.99_Z, "Spanish derived allele present in at least one copy\nat the same locus across samples in populations\n(0.99 cutoff - Z)")
  dev.off()
  
  
  
  png(filename = "vennHD_0.99.png", width = 1650, height = 2100, res = opt$res)
  plot_venn(fit_list_HD_0.99, "House derived allele present in at least one copy\nat the same locus across samples in populations\n(0.99 cutoff)")
  dev.off()
  
  png(filename = "vennHD_0.99_autosomes.png", width = 1650, height = 2100, res = opt$res)
  plot_venn(fit_list_HD_0.99_autosomes, "House derived allele present in at least one copy\nat the same locus across samples in populations\n(0.99 cutoff - autosomes)")
  dev.off()
  
  png(filename = "vennHD_0.99_Z.png", width = 1650, height = 2100, res = opt$res)
  plot_venn(fit_list_HD_0.99_Z, "House derived allele present in at least one copy\nat the same locus across samples in populations\n(0.99 cutoff - Z)")
  dev.off()


# If 10 random corsicans, re-do subsampling 10x!!####
  
  if (toupper(opt$tenRandomCorsicans) == "YES") {
    
    # redefine the read_desFile function! We need to resample again so we can't do that on the already subsetted_df
    read_desFile_resample <- function(file_path) {
      dat <- read.csv(file_path, header = TRUE, sep = "\t") # will replace the hash before Locus with X.
      dat <- dat %>%
        rename(
          "locus" = "X.Locus", # rename column to just "locus"
          "K034" = "X034" # rename column to "K034" (the K was dropped by vcftools when generating VCFs)
        )
      return(dat)
    }
    
    cutoff_0.90_sansSAHA_forResample_dat <- read_desFile_resample(opt$sansSAHA_0.90)
    
    cutoff_0.99_sansSAHA_forResample_dat <- read_desFile_resample(opt$sansSAHA_0.99)
    
    #corsica_cutoff_resample = cutoff_0.90_sansSAHA_dat %>% # doesn't matter 0.90 or 0.99, just using to pull names
      #select(starts_with("K0"))
    
    #ten_randomCorsicans_sample2 <- sample(colnames(corsica_cutoff_resample)[1:21], 10)
    ten_randomCorsicans_sample2 <- c("K071", "K026", "K089", "K010", "K090", "K029", "K073", "K034", "K039", "K085")
    
    #ten_randomCorsicans_sample3 <- sample(colnames(corsica_cutoff_resample)[1:21], 10)
    ten_randomCorsicans_sample3 <- c("K034", "K022", "K006", "K089", "K029", "K019", "K031", "K011", "K010", "K032")
    
    #ten_randomCorsicans_sample4 <- sample(colnames(corsica_cutoff_resample)[1:21], 10)
    ten_randomCorsicans_sample4 <- c("K019", "K015", "K010", "K011", "K021", "K022", "K089", "K035", "K071", "K032")
    
    #ten_randomCorsicans_sample5 <- sample(colnames(corsica_cutoff_resample)[1:21], 10)
    ten_randomCorsicans_sample5 <- c("K010", "K021", "K034", "K073", "K090", "K029", "K026", "K039", "K019", "K015")
    
    #ten_randomCorsicans_sample6 <- sample(colnames(corsica_cutoff_resample)[1:21], 10)
    ten_randomCorsicans_sample6 <- c("K089", "K032", "K083", "K022", "K019", "K021", "K011", "K034", "K015", "K010")
    
    #ten_randomCorsicans_sample7 <- sample(colnames(corsica_cutoff_resample)[1:21], 10)
    ten_randomCorsicans_sample7 <- c("K083", "K010", "K022", "K019", "K035", "K031", "K071", "K085", "K034", "K029")
    
    #ten_randomCorsicans_sample8 <- sample(colnames(corsica_cutoff_resample)[1:21], 10)
    ten_randomCorsicans_sample8 <- c("K011", "K073", "K026", "K083", "K015", "K039", "K034", "K021", "K035", "K010")
    
    #ten_randomCorsicans_sample9 <- sample(colnames(corsica_cutoff_resample)[1:21], 10)
    ten_randomCorsicans_sample9 <- c("K019", "K029", "K032", "K010", "K071", "K026", "K006", "K073", "K035", "K089")
    
    #ten_randomCorsicans_sample10 <- sample(colnames(corsica_cutoff_resample)[1:21], 10)
    ten_randomCorsicans_sample10 <- c("K026", "K011", "K019", "K039", "K015", "K035", "K032", "K006", "K096", "K073")
    
    # Function to create dataframes of 10 random corsicans & their genotypes
    corsicaRandom_resample_df <- function(dat, sampleNumber) {
      randomResample_df <- dat %>%
        select(all_of(sampleNumber)) %>%
        mutate(across(everything(), ~ gsub("SD", "D", .)),
               across(everything(), ~ gsub("HD", "D", .))) %>% # need to check D/D, D/A, and A/A counts!
        mutate_all(~ replace(., . == "A/D", "D/A")) %>% # both are heterozygotes, just clump together
        pivot_longer(cols = everything(), names_to = "sample", values_to = "genotype") %>%
        count(genotype) %>%
        pivot_wider(names_from = genotype, values_from = n, values_fill = 0)
        
    }
    
    corsicaRandom_resample1_cutoff_0.90 <- corsicaRandom_resample_df(cutoff_0.90_sansSAHA_forResample_dat, ten_randomCorsicans_sample1)
    corsicaRandom_resample2_cutoff_0.90 <- corsicaRandom_resample_df(cutoff_0.90_sansSAHA_forResample_dat, ten_randomCorsicans_sample2)
    corsicaRandom_resample3_cutoff_0.90 <- corsicaRandom_resample_df(cutoff_0.90_sansSAHA_forResample_dat, ten_randomCorsicans_sample3)
    corsicaRandom_resample4_cutoff_0.90 <- corsicaRandom_resample_df(cutoff_0.90_sansSAHA_forResample_dat, ten_randomCorsicans_sample4)
    corsicaRandom_resample5_cutoff_0.90 <- corsicaRandom_resample_df(cutoff_0.90_sansSAHA_forResample_dat, ten_randomCorsicans_sample5)
    corsicaRandom_resample6_cutoff_0.90 <- corsicaRandom_resample_df(cutoff_0.90_sansSAHA_forResample_dat, ten_randomCorsicans_sample6)
    corsicaRandom_resample7_cutoff_0.90 <- corsicaRandom_resample_df(cutoff_0.90_sansSAHA_forResample_dat, ten_randomCorsicans_sample7)
    corsicaRandom_resample8_cutoff_0.90 <- corsicaRandom_resample_df(cutoff_0.90_sansSAHA_forResample_dat, ten_randomCorsicans_sample8)
    corsicaRandom_resample9_cutoff_0.90 <- corsicaRandom_resample_df(cutoff_0.90_sansSAHA_forResample_dat, ten_randomCorsicans_sample9)
    corsicaRandom_resample10_cutoff_0.90 = corsicaRandom_resample_df(cutoff_0.90_sansSAHA_forResample_dat, ten_randomCorsicans_sample10)
    
    corsicaRandom_resample1_cutoff_0.99 <- corsicaRandom_resample_df(cutoff_0.99_sansSAHA_forResample_dat, ten_randomCorsicans_sample1)
    corsicaRandom_resample2_cutoff_0.99 <- corsicaRandom_resample_df(cutoff_0.99_sansSAHA_forResample_dat, ten_randomCorsicans_sample2)
    corsicaRandom_resample3_cutoff_0.99 <- corsicaRandom_resample_df(cutoff_0.99_sansSAHA_forResample_dat, ten_randomCorsicans_sample3)
    corsicaRandom_resample4_cutoff_0.99 <- corsicaRandom_resample_df(cutoff_0.99_sansSAHA_forResample_dat, ten_randomCorsicans_sample4)
    corsicaRandom_resample5_cutoff_0.99 <- corsicaRandom_resample_df(cutoff_0.99_sansSAHA_forResample_dat, ten_randomCorsicans_sample5)
    corsicaRandom_resample6_cutoff_0.99 <- corsicaRandom_resample_df(cutoff_0.99_sansSAHA_forResample_dat, ten_randomCorsicans_sample6)
    corsicaRandom_resample7_cutoff_0.99 <- corsicaRandom_resample_df(cutoff_0.99_sansSAHA_forResample_dat, ten_randomCorsicans_sample7)
    corsicaRandom_resample8_cutoff_0.99 <- corsicaRandom_resample_df(cutoff_0.99_sansSAHA_forResample_dat, ten_randomCorsicans_sample8)
    corsicaRandom_resample9_cutoff_0.99 <- corsicaRandom_resample_df(cutoff_0.99_sansSAHA_forResample_dat, ten_randomCorsicans_sample9)
    corsicaRandom_resample10_cutoff_0.99 = corsicaRandom_resample_df(cutoff_0.99_sansSAHA_forResample_dat, ten_randomCorsicans_sample10)
    
    # Merge the dataframes with A/A, D/D, and D/A counts for each random sample into one
    randomResample_0.90_dfs = list(corsicaRandom_resample1_cutoff_0.90, corsicaRandom_resample2_cutoff_0.90,
                                   corsicaRandom_resample3_cutoff_0.90, corsicaRandom_resample4_cutoff_0.90,
                                   corsicaRandom_resample5_cutoff_0.90, corsicaRandom_resample6_cutoff_0.90,
                                   corsicaRandom_resample7_cutoff_0.90, corsicaRandom_resample8_cutoff_0.90,
                                   corsicaRandom_resample9_cutoff_0.90, corsicaRandom_resample10_cutoff_0.90)
    
    randomResample_0.99_dfs = list(corsicaRandom_resample1_cutoff_0.99, corsicaRandom_resample2_cutoff_0.99,
                                   corsicaRandom_resample3_cutoff_0.99, corsicaRandom_resample4_cutoff_0.99,
                                   corsicaRandom_resample5_cutoff_0.99, corsicaRandom_resample6_cutoff_0.99,
                                   corsicaRandom_resample7_cutoff_0.99, corsicaRandom_resample8_cutoff_0.99,
                                   corsicaRandom_resample9_cutoff_0.99, corsicaRandom_resample10_cutoff_0.99)
    
    # Function to merge and create plot
    mergeResamples_plot <- function(randomResample_dfs_list, title) { # title is a string and should be in quotes
      
      # Merge all the resample counts dataframes into one dataframe
      corsicaRandomResample_df <- bind_rows(randomResample_dfs_list) %>%
        mutate(resampleNumber = row_number()) %>%
        relocate(resampleNumber, .before = `-`)
      
      corsicaRandomResample_long <- corsicaRandomResample_df %>%
        pivot_longer(cols = -resampleNumber, names_to = "genotype", values_to = "count")
      
      # Create a long format dataframe for the mean, SE, and 95% confidence intervals
      mean_counts <- corsicaRandomResample_long %>%
        group_by(genotype) %>%
        summarize(
          mean_count = mean(count), # mean of column (mean genotype count across samples)
          se = sd(count) / sqrt(n()),  # sd(cosricaRandomResample_cutoff_0.90$genotypeColumn)
          lower_CI = mean_count - qt(0.975, df = n() - 1) * se, # Rbloggers, Calculate Confidence Intervals in R
          upper_CI = mean_count + qt(0.975, df = n() - 1) * se # Confidence Interval = (point estimate)+/-(critical value)*(standard error)
        )
      
      # Create the bar plot
      plot <- ggplot(corsicaRandomResample_long, aes(x = genotype, y = count, fill = as.factor(resampleNumber))) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_minimal() +
        labs(x = "Genotype",
             y = "Genotype Count",
             fill = "Resample number") +
        ggtitle(title) +
        scale_y_continuous(labels = scales::comma_format()) +
        geom_line(data = mean_counts, aes(x = genotype, y = mean_count, group = 1), 
                  color = "black", linetype = "dashed", size = 1, inherit.aes = FALSE) +
        geom_point(data = mean_counts, aes(x = genotype, y = mean_count), 
                   color = "black", size = 3, inherit.aes = FALSE) +
        geom_errorbar(data = mean_counts, aes(x = genotype, ymin = lower_CI, ymax = upper_CI), 
                      width = 0.2, color = "black", inherit.aes = FALSE)
      
      return(plot)
      
    }
    
    plot_0.90 <- mergeResamples_plot(randomResample_0.90_dfs, "Corsica genotype counts per resampling event\n(10 random Corsican samples selected each time)\n0.90 cutoff")
    ggsave(plot_0.90, filename = "corsicaResample_0.90.png", width = 6, height = 4, dpi = opt$res, bg = "white") # size is in pixels, not inches
    
    plot_0.99 <- mergeResamples_plot(randomResample_0.99_dfs, "Corsica genotype counts per resampling event\n(10 random Corsican samples selected each time)\n0.99 cutoff")
    ggsave(plot_0.99, filename = "corsicaResample_0.99.png", width = 5.7, height = 4, dpi = opt$res, bg = "white")
    
    # dev.off not working, use ggsave instead
    #png(filename = "corsicaResample_0.90.png", width = opt$width, height = opt$height, res = opt$res)
    #mergeResamples_plot(randomResample_0.90_dfs, "Corsica genotype counts per resampling event\n(10 random Corsican samples selected each time)\n0.90 cutoff")
    #dev.off()
    
    #png(filename = "corsicaResample_0.99.png", width = opt$width, height = opt$height, res = opt$res)
    #mergeResamples_plot(randomResample_0.99_dfs, "Corsica genotype counts per resampling event\n(10 random Corsican samples selected each time)\n0.99 cutoff")
    #dev.off()
    
  } 