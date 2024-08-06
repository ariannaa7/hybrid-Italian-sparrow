#!/usr/bin/env Rscript

# R script with output for all final tables & figures

# Example usage: 

# Rscript finalPlots.R --SAHA_0.99 designationFiles/cutoff_0.99_SAHA.tsv \
# --tenRandomCorsicans NO --lpData YES

# Install & load necessary packages ####

# List of packages needed
packages <- c("optparse", # v 1.7.5, allow use of command line flags
              "tidyverse", # v 2.0.0, data frame manipulation & plotting
              "ggpubr", # v. 0.6.0, arrange multiple plots togehter, chi-sq observed-expected plot
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
  
  make_option(c("--SAHA_0.99"), type = "character", default = NULL, 
              help = "Input designation file for 0.99 cutoff, SAHA", metavar = "file_path"),
  
  make_option(c("--tenRandomCorsicans"), type = "character", default = "NO", 
              help = "write 'yes' or 'no'. The analysis should proceed with the 10 random corsicans = 'yes'. ", metavar = "character"),
  
  make_option(c("--lpData"), type = "character", default = "NO", 
              help = "write 'yes' or 'no'. You've provided linkage pruned data = 'yes'. ", metavar = "character"),
  
  make_option(c("--res"), type = "integer", default = 300, 
              help = "Resolution of output images in DPI [default= %default]", metavar = "integer")
  
)

# Set the flag list
opt_parser <- OptionParser(option_list = flag_list)
opt <- parse_args(opt_parser) # can refrence user provided arguments using opt$res for example!

# Import Data ####

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
cutoff_0.99_SAHA_dat <- read_desFile(opt$SAHA_0.99) # includes autosomes & Z, with or without corsicans depending on user choice
cutoff_0.99_SAHA_dat_autosomes <- desFile_autosomes(cutoff_0.99_SAHA_dat) # just autosomes
cutoff_0.99_SAHA_dat_Z <- desFile_Z(cutoff_0.99_SAHA_dat) # just Z

# Italian heterozygote/homozygote frequency ####

# Create a function that will make dataframe with all the italian genotypes in the designation file grouped by population
genotypeFreqPlot <- function(dat) {
  
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
  
  italGenotypesCounts_df <- ital_genotypes_df %>%
    pivot_longer(cols = everything(), names_to = "population", values_to = "genotype") %>% # transform to long format
    group_by(population, genotype) %>% # group by population and genotype
    summarize(count = n()) %>% # then count occurrence of each genotype for each population
    pivot_wider(names_from = genotype, values_from = count, values_fill = list(count = 0)) %>% # change back to wide format to easily remove NA column
    select(-"NA") # remove the NA column, only exists because more K0 samples than all others
  
  zygosity_counts_SAHA_df <- italGenotypesCounts_df %>%
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
  
  long_format <- zygosity_freq_SAHA_df %>%
    pivot_longer(cols = starts_with("freq"), names_to = "genotype", values_to = "frequency")
  
  plot <- ggplot(long_format, aes(x = population, y = frequency, fill = genotype)) +
    geom_bar(stat = "identity") +
    labs(x = "Italian population", y = "Genotype frequency", fill = "") + # don't want a legend title
    scale_fill_manual(values = c("freq_missing_genotype" = "lawngreen", "freq_homozygous_houseAncestral" = "#5c2be2", "freq_homozygous_spanishAncestral" = "magenta2", "freq_homozygous_houseDerived" = "royalblue", "freq_homozygous_spanishDerived" = "orangered3", "freq_heterozygous_houseDerived_spanishAncestral" = "skyblue1", "freq_heterozygous_spanishDerived_houseAncestral" = "sienna2"), 
                      labels = c("freq_missing_genotype" = "Missing genotype", "freq_homozygous_houseAncestral" = "Homozygous: House Ancestral", "freq_homozygous_spanishAncestral" = "Homozygous: Spanish Ancestral", "freq_homozygous_houseDerived" = "Homozgyous: House Derived", "freq_homozygous_spanishDerived" = "Homozygous: Spanish Derived", "freq_heterozygous_houseDerived_spanishAncestral" = "Heterozygous: House Derived/Spanish Ancestral", "freq_heterozygous_spanishDerived_houseAncestral" = "Heterozygous: Spanish Derived/House Ancestral")) +
    scale_x_discrete(labels = c("C0" = "Crete", "K0" = "Corsica", "M0" = "Malta", "S0" = "Sicily" )) +
    #ggtitle(title) +
    theme_minimal()

  return(plot)
}

# Plot the charts using the functions
png(filename = "genotypeFreqPerPop_0.99_SAHA.png", width = 2000, height = 1300, res = opt$res)
genotypeFreqPlot(cutoff_0.99_SAHA_dat) # Genotype frequency per population at 0.99 cutoff
dev.off()

png(filename = "genotypeFreqPerPop_0.99_SAHA_autosomes.png", width = 2000, height = 1300, res = opt$res)
genotypeFreqPlot(cutoff_0.99_SAHA_dat_autosomes) # Genotype frequency per population at 0.99 cutoff - autosomes
dev.off()

png(filename = "genotypeFreqPerPop_0.99_SAHA_Z.png", width = 2000, height = 1300, res = opt$res)
genotypeFreqPlot(cutoff_0.99_SAHA_dat_Z) # Genotype frequency per population at 0.99 cutoff - Z
dev.off()


# Venn diagrams! ####

# For each population, create a list with locus and allele (split up genotype) for every sample
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

corsica_vennList_cutoff_0.99 <- vennList(cutoff_0.99_SAHA_dat, popPrefix = "K0")

crete_vennList_cutoff_0.99 <- vennList(cutoff_0.99_SAHA_dat, popPrefix = "C0")

malta_vennList_cutoff_0.99 <- vennList(cutoff_0.99_SAHA_dat, popPrefix = "M0")

sicily_vennList_cutoff_0.99 <- vennList(cutoff_0.99_SAHA_dat, popPrefix = "S0")

# Use eulerr package to create a fit list

# Look at any occurence of SD at the locus (can be heterozygote/homozygote)
fit_list_SD_0.99 <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.99$locus[corsica_vennList_cutoff_0.99$allele == "SD"]),
                                      Crete = unique(crete_vennList_cutoff_0.99$locus[crete_vennList_cutoff_0.99$allele == "SD"]),
                                      Malta = unique(malta_vennList_cutoff_0.99$locus[malta_vennList_cutoff_0.99$allele == "SD"]),
                                      Sicily = unique(sicily_vennList_cutoff_0.99$locus[sicily_vennList_cutoff_0.99$allele == "SD"]))
)

# Look at any occurence of HD at the locus (can be heterozygote/homozygote)
fit_list_HD_0.99 <- eulerr::venn(list(Corsica = unique(corsica_vennList_cutoff_0.99$locus[corsica_vennList_cutoff_0.99$allele == "HD"]),
                                      Crete = unique(crete_vennList_cutoff_0.99$locus[crete_vennList_cutoff_0.99$allele == "HD"]),
                                      Malta = unique(malta_vennList_cutoff_0.99$locus[malta_vennList_cutoff_0.99$allele == "HD"]),
                                      Sicily = unique(sicily_vennList_cutoff_0.99$locus[sicily_vennList_cutoff_0.99$allele == "HD"]))
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

# Spanish derived allele present in at least one copy at the same locus across samples in populations (0.99 cutoff)
png(filename = "vennSD_0.99.png", width = 1700, height = 2100, res = opt$res)
plot_venn(fit_list_SD_0.99, "")
dev.off()

# House derived allele present in at least one copy at the same locus across samples in populations (0.99 cutoff)
png(filename = "vennHD_0.99.png", width = 1650, height = 2100, res = opt$res)
plot_venn(fit_list_HD_0.99, "")
dev.off()


# Counts - where in genome, autosome & Z ####

VCF_intersection_stats = data.frame( 
  chrom_type = c("Autosomes", "Z"),
  exon_sites_count = c(2271,6905),
  gene_sites_count = c(17473,64307),
  intergenic_sites_count = c(7600,48031),
  nearGene_sites_count = c(3525,20467) # within 10,000 bps of gene
)

# Long format
VCF_intersection_stats_long <- pivot_longer(VCF_intersection_stats, cols = -chrom_type, names_to = "count_type", values_to = "count")

# Plot
VCF_intersection_stats_plot <-  ggplot(VCF_intersection_stats_long, aes(x = as.factor(chrom_type), y = count, fill = count_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "", fill = "Feature Type") + # y = Number of designated sites falling within this region
  scale_y_continuous(labels = scales::comma_format()) +
  scale_fill_manual(values = c("exon_sites_count" = "seagreen", "gene_sites_count" = "sienna", "intergenic_sites_count" = "slateblue", "nearGene_sites_count" = "palevioletred"),
                    labels = c("exon_sites_count" = "Exons", "gene_sites_count" = "Genes", "intergenic_sites_count" = "Intergenic", "nearGene_sites_count" = "Near gene")) +
  theme_minimal() +
  theme(legend.title = element_text(hjust = 0.5), # center the legend title
        axis.text.x = element_text(size = 11)) # increase x-categorical label sizes

# Save output
png(filename = "SiteRegions_chrom_noPrune.png", width = 2000, height = 2000, res = opt$res)

VCF_intersection_stats_plot

dev.off()

# Test #1 - One sample proportion test, test for tendency to generally purge derived sites (test derived against expectation of 50%) ####

# Function that will create a dataframe with D, A, total allele counts for each Italian population
count_Italian_DA <- function(dat, population_prefix) { # population_prefix should be entered as a string in ""
  DA_count <- dat %>% # dat being the designation dataframe
    select(starts_with(population_prefix)) %>% # only pull columns from the designation df which start with the provided population prefix
    pivot_longer(everything(), names_to = "sample", values_to = "genotype") %>% # change to long format 
    group_by(sample) %>% # group by sample, (e.g. corsica 0.99 - 2,885,631 rows at this point bc 21 samples * 137,411 loci)
    filter(genotype != "-") %>% # remove missing genotypes, (e.g. corsica 0.99 - 2,878,328 rows at this point)
    summarize(D_count = sum(str_count(genotype, "D")), # will count instances of "D" present in each genotype etc.
              A_count = sum(str_count(genotype, "A")))
  
  # Create a dataframe with count totals for the whole population, not just for samples within the population
  DA_count_totals <- DA_count %>%
    summarize(
      D_total = sum(D_count), # (e.g. corsica 0.99 - 2,745,426 alleles)
      A_total = sum(A_count) # (e.g. corsica 0.99 - 3,011,230 alleles)
    ) 
  
  # Add on columns for House total, spanish total, and total total
  DA_count_totals <- DA_count_totals %>%
    mutate(totaAlleles = D_total + A_total) # e.g. corsica 0.99 sum = 5,756,656 = (2,878,328 rows * 2 alleles per row)
  
  return(DA_count_totals)
}

# Run the function
corsica_DA_Counts_0.99 <- count_Italian_DA(cutoff_0.99_SAHA_dat, population_prefix = "K0")
corsica_DA_Counts_0.99_autosomes <- count_Italian_DA(cutoff_0.99_SAHA_dat_autosomes, population_prefix = "K0")
corsica_DA_Counts_0.99_Z <- count_Italian_DA(cutoff_0.99_SAHA_dat_Z, population_prefix = "K0")

crete_DA_Counts_0.99 <- count_Italian_DA(cutoff_0.99_SAHA_dat, population_prefix = "C0")
crete_DA_Counts_0.99_autosomes <- count_Italian_DA(cutoff_0.99_SAHA_dat_autosomes, population_prefix = "C0")
crete_DA_Counts_0.99_Z <- count_Italian_DA(cutoff_0.99_SAHA_dat_Z, population_prefix = "C0")

malta_DA_Counts_0.99 <- count_Italian_DA(cutoff_0.99_SAHA_dat, population_prefix = "M0")
malta_DA_Counts_0.99_autosomes <- count_Italian_DA(cutoff_0.99_SAHA_dat_autosomes, population_prefix = "M0")
malta_DA_Counts_0.99_Z <- count_Italian_DA(cutoff_0.99_SAHA_dat_Z, population_prefix = "M0")

sicily_DA_Counts_0.99 <- count_Italian_DA(cutoff_0.99_SAHA_dat, population_prefix = "S0")
sicily_DA_Counts_0.99_autosomes <- count_Italian_DA(cutoff_0.99_SAHA_dat_autosomes, population_prefix = "S0")
sicily_DA_Counts_0.99_Z <- count_Italian_DA(cutoff_0.99_SAHA_dat_Z, population_prefix = "S0")


# Create a prop test function
prop_test <- function(population_DA_Counts_cutoff) {
  prop <- prop.test(x = population_DA_Counts_cutoff$D_total, # number of derived alleles
                    n = population_DA_Counts_cutoff$totaAlleles, # number of alleles (D + A)
                    p = 0.5, # probability we are testing against, no purging means 50% derived 50% ancestral
                    correct = FALSE,
                    alternative = "two.sided" # null (H0): p = 0.50, alternate (Ha): p < 0.50
  )
  return(prop)
}

# Run the function and save the output!
sink("corsica_0.99_prop.test.out")
corsica_prop_0.99 <- prop_test(corsica_DA_Counts_0.99) # run test
print(corsica_prop_0.99) # save output
sink()

sink("corsica_0.99_prop.test_autosomes.out")
corsica_prop_0.99_autosomes <- prop_test(corsica_DA_Counts_0.99_autosomes)
print(corsica_prop_0.99_autosomes)
sink()

sink("corsica_0.99_prop.test_Z.out")
corsica_prop_0.99_Z <- prop_test(corsica_DA_Counts_0.99_Z)
print(corsica_prop_0.99_Z)
sink()


sink("crete_0.99_prop.test.out")
crete_prop_0.99 <- prop_test(crete_DA_Counts_0.99)
print(crete_prop_0.99)
sink()

sink("crete_0.99_prop.test_autosomes.out")
crete_prop_0.99_autosomes <- prop_test(crete_DA_Counts_0.99_autosomes)
print(crete_prop_0.99_autosomes)
sink()

sink("crete_0.99_prop.test_Z.out")
crete_prop_0.99_Z <- prop_test(crete_DA_Counts_0.99_Z)
print(crete_prop_0.99_Z)
sink()


sink("malta_0.99_prop.test.out")
malta_prop_0.99 <- prop_test(malta_DA_Counts_0.99)
print(malta_prop_0.99)
sink()

sink("malta_0.99_prop.test_autosomes.out")
malta_prop_0.99_autosomes <- prop_test(malta_DA_Counts_0.99_autosomes)
print(malta_prop_0.99_autosomes)
sink()

sink("malta_0.99_prop.test_Z.out")
malta_prop_0.99_Z <- prop_test(malta_DA_Counts_0.99_Z)
print(malta_prop_0.99_Z)
sink()


sink("sicily_0.99_prop.test.out")
sicily_prop_0.99 <- prop_test(sicily_DA_Counts_0.99)
print(sicily_prop_0.99)
sink()

sink("sicily_0.99_prop.test_autosomes.out")
sicily_prop_0.99_autosomes <- prop_test(sicily_DA_Counts_0.99_autosomes)
print(sicily_prop_0.99_autosomes)
sink()

sink("sicily_0.99_prop.test_Z.out")
sicily_prop_0.99_Z <- prop_test(sicily_DA_Counts_0.99_Z)
print(sicily_prop_0.99_Z)
sink()


# Test #2 - Chisq (ind) HA/HD/SA/SD test if proportion derived differs between house and Spanish ####
  # Create allele counts dataframe - for manual calc ####
  
  # Function that will create a dataframe with HD, SD, HA, SA, house, spanish, and total allele counts for each Italian population
  count_Italian_alleles <- function(dat, population_prefix) { # population_prefix should be entered as a string in ""
    allele_count <- dat %>% # dat being the designation dataframe
      select(starts_with(population_prefix)) %>% # only pull columns from the designation df which start with the provided population prefix
      pivot_longer(everything(), names_to = "sample", values_to = "genotype") %>% # change to long format 
      group_by(sample) %>% # group by sample, (e.g. corsica 0.99 - 2,885,631 rows at this point bc 21 samples * 137,411 loci)
      filter(genotype != "-") %>% # remove missing genotypes, (e.g. corsica 0.99 - 2,878,328 rows at this point)
      summarize(HD_count = sum(str_count(genotype, "HD")), # will count instances of "HD" present in each genotype etc.
                SD_count = sum(str_count(genotype, "SD")),
                HA_count = sum(str_count(genotype, "HA")),
                SA_count = sum(str_count(genotype, "SA")))
    
    # Create a dataframe with count totals for the whole population, not just for samples within the population
    allele_count_totals <- allele_count %>%
      summarize(
        HD_total = sum(HD_count), # (e.g. corsica 0.99 - 1,710,113 alleles)
        SD_total = sum(SD_count), # (e.g. corsica 0.99 - 1,035,313 alleles)
        HA_total = sum(HA_count), # (e.g. corsica 0.99 - 2,200,151 alleles)
        SA_total = sum(SA_count) # (e.g. corsica 0.99 - 811,079 alleles)
      ) 
    
    # Add on columns for House total, spanish total, and total total
    allele_count_totals <- allele_count_totals %>%
      mutate(totalHouseAlleles = HD_total + HA_total) %>%
      mutate(totalSpanishAlleles = SD_total + SA_total) %>%
      mutate(totaAlleles = HD_total + SD_total + HA_total + SA_total) # e.g. corsica 0.99 sum = 5,756,656 = (2,878,328 rows * 2 alleles per row)
    
    return(allele_count_totals)
  }

  # Run the function for each population
  corsica_alleleCounts_0.99 <- count_Italian_alleles(cutoff_0.99_SAHA_dat, population_prefix = "K0")
  corsica_alleleCounts_0.99_autosomes <- count_Italian_alleles(cutoff_0.99_SAHA_dat_autosomes, population_prefix = "K0")
  corsica_alleleCounts_0.99_Z <- count_Italian_alleles(cutoff_0.99_SAHA_dat_Z, population_prefix = "K0")
  
  crete_alleleCounts_0.99 <- count_Italian_alleles(cutoff_0.99_SAHA_dat, population_prefix = "C0")
  crete_alleleCounts_0.99_autosomes <- count_Italian_alleles(cutoff_0.99_SAHA_dat_autosomes, population_prefix = "C0")
  crete_alleleCounts_0.99_Z <- count_Italian_alleles(cutoff_0.99_SAHA_dat_Z, population_prefix = "C0")
  
  malta_alleleCounts_0.99 <- count_Italian_alleles(cutoff_0.99_SAHA_dat, population_prefix = "M0")
  malta_alleleCounts_0.99_autosomes <- count_Italian_alleles(cutoff_0.99_SAHA_dat_autosomes, population_prefix = "M0")
  malta_alleleCounts_0.99_Z <- count_Italian_alleles(cutoff_0.99_SAHA_dat_Z, population_prefix = "M0")
  
  sicily_alleleCounts_0.99 <- count_Italian_alleles(cutoff_0.99_SAHA_dat, population_prefix = "S0")
  sicily_alleleCounts_0.99_autosomes <- count_Italian_alleles(cutoff_0.99_SAHA_dat_autosomes, population_prefix = "S0")
  sicily_alleleCounts_0.99_Z <- count_Italian_alleles(cutoff_0.99_SAHA_dat_Z, population_prefix = "S0")


  # Create chi-sq dataframe - for manual calc ####
  
  # Function that creates chi-sq dataframe and computes chi-sq value manually
  chi_squared_test2_calc <- function(alleleCounts_df) { # population_prefix should be entered as a string in ""
    
    chi_squared_df <- data.frame(
      combination = c("House_derived", "House_ancestral", "Spanish_derived", "Spanish_ancestral")
    )
    
    # e.g. ((# House derived alleles + # Spanish derived alleles)/ # allles total)
    propDerived_value <- ((alleleCounts_df$HD_total + alleleCounts_df$SD_total)/alleleCounts_df$totaAlleles)
    
    
    propAncestral_value <- (1-propDerived_value)
    
    # add a column to df that shows the expected proportion
    chi_squared_df <- chi_squared_df %>% # proportion ancestral is 1-propDerived_value
      bind_cols(expected_prop = c(propDerived_value, propAncestral_value, propDerived_value, propAncestral_value))
    
    # add a column to df that shows the expected number
    chi_squared_df <- chi_squared_df %>%
      bind_cols(expected_number = c((propDerived_value * alleleCounts_df$totalHouseAlleles), (propAncestral_value * alleleCounts_df$totalHouseAlleles),
                                    (propDerived_value * alleleCounts_df$totalSpanishAlleles), (propAncestral_value * alleleCounts_df$totalSpanishAlleles)))
    
    # add a column to df that shows the observed number
    chi_squared_df <- chi_squared_df %>%
      bind_cols(observed_number = c(alleleCounts_df$HD_total, alleleCounts_df$HA_total, alleleCounts_df$SD_total, alleleCounts_df$SA_total))
    
    # add a column to df that shows the result of (observed - expected)^2 / expected
    chi_squared_df <- chi_squared_df %>%
      bind_cols(chiSquared = (data.frame(chiSquared = (.$observed_number - .$expected_number)^2 / .$expected_number)))
    
    chi_squared_row <- c(NA, NA, NA, NA, sum(chi_squared_df$chiSquared)) # sum the last column to get chi squared value
    
    # Append the last row to the dataframe
    chi_squared_df <- rbind(chi_squared_df, chi_squared_row)
    
    return(chi_squared_df)
  }
  
  # Run the function for each population & save output
  corsica_chiSquared2_df_0.99 <- chi_squared_test2_calc(corsica_alleleCounts_0.99)
  write.csv(corsica_chiSquared2_df_0.99, "corsica_chiSquared2_df_0.99.csv")
  
  corsica_chiSquared2_df_0.99_autosomes <- chi_squared_test2_calc(corsica_alleleCounts_0.99_autosomes)
  write.csv(corsica_chiSquared2_df_0.99_autosomes, "corsica_chiSquared2_df_0.99_autosomes.csv")
  
  corsica_chiSquared2_df_0.99_Z <- chi_squared_test2_calc(corsica_alleleCounts_0.99_Z)
  write.csv(corsica_chiSquared2_df_0.99_Z, "corsica_chiSquared2_df_0.99_Z.csv")
  
  
  crete_chiSquared2_df_0.99 <- chi_squared_test2_calc(crete_alleleCounts_0.99)
  write.csv(crete_chiSquared2_df_0.99, "crete_chiSquared2_df_0.99.csv")
  
  crete_chiSquared2_df_0.99_autosomes <- chi_squared_test2_calc(crete_alleleCounts_0.99_autosomes)
  write.csv(crete_chiSquared2_df_0.99_autosomes, "crete_chiSquared2_df_0.99_autosomes.csv")
  
  crete_chiSquared2_df_0.99_Z <- chi_squared_test2_calc(crete_alleleCounts_0.99_Z)
  write.csv(crete_chiSquared2_df_0.99_Z, "crete_chiSquared2_df_0.99_Z.csv")
  
  
  malta_chiSquared2_df_0.99 <- chi_squared_test2_calc(malta_alleleCounts_0.99)
  write.csv(malta_chiSquared2_df_0.99, "malta_chiSquared2_df_0.99.csv")
  
  malta_chiSquared2_df_0.99_autosomes <- chi_squared_test2_calc(malta_alleleCounts_0.99_autosomes)
  write.csv(malta_chiSquared2_df_0.99_autosomes, "malta_chiSquared2_df_0.99_autosomes.csv")
  
  malta_chiSquared2_df_0.99_Z <- chi_squared_test2_calc(malta_alleleCounts_0.99_Z)
  write.csv(malta_chiSquared2_df_0.99_Z, "malta_chiSquared2_df_0.99_Z.csv")
  
  
  sicily_chiSquared2_df_0.99 <- chi_squared_test2_calc(sicily_alleleCounts_0.99)
  write.csv(sicily_chiSquared2_df_0.99, "sicily_chiSquared2_df_0.99.csv")
  
  sicily_chiSquared2_df_0.99_autosomes <- chi_squared_test2_calc(sicily_alleleCounts_0.99_autosomes)
  write.csv(sicily_chiSquared2_df_0.99_autosomes, "sicily_chiSquared2_df_0.99_autosomes.csv")
  
  sicily_chiSquared2_df_0.99_Z <- chi_squared_test2_calc(sicily_alleleCounts_0.99_Z)
  write.csv(sicily_chiSquared2_df_0.99_Z, "sicily_chiSquared2_df_0.99_Z.csv")
  
  # Create contingency dataframe - for built in calc ####
  
  # Function to create contignency "tables" (actually using df format), Reference: http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
  contigency_test2_df <-function(population_alleleCounts_cutoff) {
    
    count_matrix <- matrix(c(population_alleleCounts_cutoff$HA_total, population_alleleCounts_cutoff$HD_total, population_alleleCounts_cutoff$SA_total, population_alleleCounts_cutoff$SD_total), nrow = 2, byrow = TRUE,
                           dimnames = list(c("house", "spanish"),
                                           c("ancestral", "derived")))
    return(count_matrix)
  }
  
  # Run the function to create contingency tables
  corsica_contigency2_0.99 <- contigency_test2_df(corsica_alleleCounts_0.99)
  corsica_contigency2_0.99_autosomes <- contigency_test2_df(corsica_alleleCounts_0.99_autosomes)
  corsica_contigency2_0.99_Z <- contigency_test2_df(corsica_alleleCounts_0.99_Z)
  
  crete_contigency2_0.99 <- contigency_test2_df(crete_alleleCounts_0.99)
  crete_contigency2_0.99_autosomes <- contigency_test2_df(crete_alleleCounts_0.99_autosomes)
  crete_contigency2_0.99_Z <- contigency_test2_df(crete_alleleCounts_0.99_Z)
  
  malta_contigency2_0.99 <- contigency_test2_df(malta_alleleCounts_0.99)
  malta_contigency2_0.99_autosomes <- contigency_test2_df(malta_alleleCounts_0.99_autosomes)
  malta_contigency2_0.99_Z <- contigency_test2_df(malta_alleleCounts_0.99_Z)
  
  sicily_contigency2_0.99 <- contigency_test2_df(sicily_alleleCounts_0.99)
  sicily_contigency2_0.99_autosomes <- contigency_test2_df(sicily_alleleCounts_0.99_autosomes)
  sicily_contigency2_0.99_Z <- contigency_test2_df(sicily_alleleCounts_0.99_Z)
  
  # Run chisq.test & save output - for built-in calc ####
  
  # This function uses contingency tables to calc chi-sq, calculates expected by doing (row_total * column_total)/total, same expected as if done manually (achieves same thing)
  sink("corsica_chisq2_0.99.out")
  corsica_chisq <- chisq.test(corsica_contigency2_0.99)
  print(corsica_chisq)
  sink()
  
  sink("corsica_chisq2_auto_0.99.out")
  corsica_chisq_auto <- chisq.test(corsica_contigency2_0.99_autosomes)
  print(corsica_chisq_auto)
  sink()
  
  sink("corsica_chisq2_Z_0.99.out")
  corsica_chisq_Z <- chisq.test(corsica_contigency2_0.99_Z)
  print(corsica_chisq_Z)
  sink()
  
  
  sink("crete_chisq2_0.99.out")
  crete_chisq <- chisq.test(crete_contigency2_0.99)
  print(crete_chisq)
  sink()
  
  sink("crete_chisq2_auto_0.99.out")
  crete_chisq_auto <- chisq.test(crete_contigency2_0.99_autosomes)
  print(crete_chisq_auto)
  sink()

  sink("crete_chisq2_Z_0.99.out")
  crete_chisq_Z <- chisq.test(crete_contigency2_0.99_Z)
  print(crete_chisq_Z)
  sink()
  
  
  sink("malta_chisq2_0.99.out")
  malta_chisq <- chisq.test(malta_contigency2_0.99)
  print(malta_chisq)
  sink()
  
  sink("malta_chisq2_auto_0.99.out")
  malta_chisq_auto <- chisq.test(malta_contigency2_0.99_autosomes)
  print(malta_chisq_auto)
  sink()
  
  sink("malta_chisq2_Z_0.99.out")
  malta_chisq_Z <- chisq.test(malta_contigency2_0.99_Z)
  print(malta_chisq_Z)
  sink()
  

  sink("sicily_chisq2_0.99.out")
  sicily_chisq <- chisq.test(sicily_contigency2_0.99)
  print(sicily_chisq)
  sink()
  
  sink("sicily_chisq2_auto_0.99.out")
  sicily_chisq_auto <- chisq.test(sicily_contigency2_0.99_autosomes)
  print(sicily_chisq_auto)
  sink()
  
  sink("malta_chisq2_Z_0.99.out")
  sicily_chisq_Z <- chisq.test(sicily_contigency2_0.99_Z)
  print(sicily_chisq_Z)
  sink()
  
# Test #3 - Chisq (gof) sep by parentage, post-hoc, tests if purging is species specific ####
  # Create allele counts dataframes, sep by parent - for manual calc ####
  
  # Function to pull allele counts for Spanish df
  count_Italian_alleles_spanish <- function(allele_count_df) { # population_prefix should be entered as a string in ""
    
    spanish_df <- allele_count_df %>%
      select(SD_total, SA_total, totalSpanishAlleles)
    
    return(spanish_df)
  }
  
  # Function to pull allele counts for house df
  count_Italian_alleles_house <- function(allele_count_df) { # population_prefix should be entered as a string in ""
    
    house_df <- allele_count_df %>%
      select(HD_total, HA_total, totalHouseAlleles)
    
    return(house_df)
  }
  
  # Run the functions!
  corsica_spanish_alleleCounts_0.99_autosomes <- count_Italian_alleles_spanish(corsica_alleleCounts_0.99_autosomes)
  corsica_spanish_alleleCounts_0.99_Z <- count_Italian_alleles_spanish(corsica_alleleCounts_0.99_Z)
  corsica_house_alleleCounts_0.99_autosomes <- count_Italian_alleles_house(corsica_alleleCounts_0.99_autosomes)
  corsica_house_alleleCounts_0.99_Z <- count_Italian_alleles_house(corsica_alleleCounts_0.99_Z)
  
  crete_spanish_alleleCounts_0.99_autosomes <- count_Italian_alleles_spanish(crete_alleleCounts_0.99_autosomes)
  crete_spanish_alleleCounts_0.99_Z <- count_Italian_alleles_spanish(crete_alleleCounts_0.99_Z)
  crete_house_alleleCounts_0.99_autosomes <- count_Italian_alleles_house(crete_alleleCounts_0.99_autosomes)
  crete_house_alleleCounts_0.99_Z <- count_Italian_alleles_house(crete_alleleCounts_0.99_Z)
  
  malta_spanish_alleleCounts_0.99_autosomes <- count_Italian_alleles_spanish(malta_alleleCounts_0.99_autosomes)
  malta_spanish_alleleCounts_0.99_Z <- count_Italian_alleles_spanish(malta_alleleCounts_0.99_Z)
  malta_house_alleleCounts_0.99_autosomes <- count_Italian_alleles_house(malta_alleleCounts_0.99_autosomes)
  malta_house_alleleCounts_0.99_Z <- count_Italian_alleles_house(malta_alleleCounts_0.99_Z)
  
  sicily_spanish_alleleCounts_0.99_autosomes <- count_Italian_alleles_spanish(sicily_alleleCounts_0.99_autosomes)
  sicily_spanish_alleleCounts_0.99_Z <- count_Italian_alleles_spanish(sicily_alleleCounts_0.99_Z)
  sicily_house_alleleCounts_0.99_autosomes <- count_Italian_alleles_house(sicily_alleleCounts_0.99_autosomes)
  sicily_house_alleleCounts_0.99_Z <- count_Italian_alleles_house(sicily_alleleCounts_0.99_Z)
  
  # Create chi-sq dataframes, sep by parent - for manual calc ####
  chi_squared_test3_spanish_calc <- function(alleleCounts_df) { # population_prefix should be entered as a string in ""
    
    chi_squared_df <- data.frame(
      combination = c("Spanish_derived", "Spanish_ancestral")
    )
    
    # e.g. ((# House derived alleles + # Spanish derived alleles)/ # allles total)
    propDerived_value <- (0.50) # only two possibilites, so we can't calc a prop based on observed, go 50/50
    
    propAncestral_value <- (0.50)
    
    # add a column to df that shows the expected proportion
    chi_squared_df <- chi_squared_df %>% # proportion ancestral is 1-propDerived_value
      bind_cols(expected_prop = c(propDerived_value, propAncestral_value))
    
    # add a column to df that shows the expected number
    chi_squared_df <- chi_squared_df %>%
      bind_cols(expected_number = c((propDerived_value * alleleCounts_df$totalSpanishAlleles), (propAncestral_value * alleleCounts_df$totalSpanishAlleles)))
    
    # add a column to df that shows the observed number
    chi_squared_df <- chi_squared_df %>%
      bind_cols(observed_number = c(alleleCounts_df$SD_total, alleleCounts_df$SA_total))
    
    # add a column to df that shows the result of (observed - expected)^2 / expected
    chi_squared_df <- chi_squared_df %>%
      bind_cols(chiSquared = (data.frame(chiSquared = (.$observed_number - .$expected_number)^2 / .$expected_number)))
    
    chi_squared_row <- c(NA, NA, NA, NA, sum(chi_squared_df$chiSquared)) # sum the last column to get chi squared value
    
    # Append the last row to the dataframe
    chi_squared_df <- rbind(chi_squared_df, chi_squared_row)
    
    return(chi_squared_df)
  }
  
  chi_squared_test3_house_calc <- function(alleleCounts_df) { # population_prefix should be entered as a string in ""
    
    chi_squared_df <- data.frame(
      combination = c("House_derived", "House_ancestral")
    )
    
    propDerived_value <- (0.50) # only two possibilites, so we can't calc a prop based on observed, go 50/50
    
    propAncestral_value <- (0.50)
    
    # add a column to df that shows the expected proportion
    chi_squared_df <- chi_squared_df %>% # proportion ancestral is 1-propDerived_value
      bind_cols(expected_prop = c(propDerived_value, propAncestral_value))
    
    # add a column to df that shows the expected number
    chi_squared_df <- chi_squared_df %>%
      bind_cols(expected_number = c((propDerived_value * alleleCounts_df$totalHouseAlleles), (propAncestral_value * alleleCounts_df$totalHouseAlleles)))
    
    # add a column to df that shows the observed number
    chi_squared_df <- chi_squared_df %>%
      bind_cols(observed_number = c(alleleCounts_df$HD_total, alleleCounts_df$HA_total))
    
    # add a column to df that shows the result of (observed - expected)^2 / expected
    chi_squared_df <- chi_squared_df %>%
      bind_cols(chiSquared = (data.frame(chiSquared = (.$observed_number - .$expected_number)^2 / .$expected_number)))
    
    chi_squared_row <- c(NA, NA, NA, NA, sum(chi_squared_df$chiSquared)) # sum the last column to get chi squared value
    
    # Append the last row to the dataframe
    chi_squared_df <- rbind(chi_squared_df, chi_squared_row)
    
    return(chi_squared_df)
  }
  
  corsica_spanish_chiSquared3_df_0.99_autosomes <- chi_squared_test3_spanish_calc(corsica_spanish_alleleCounts_0.99_autosomes)
  write.csv(corsica_spanish_chiSquared3_df_0.99_autosomes, "corsica_spanish_chiSquared3_df_0.99_autosomes.csv")
  
  corsica_spanish_chiSquared3_df_0.99_Z <- chi_squared_test3_spanish_calc(corsica_spanish_alleleCounts_0.99_Z)
  write.csv(corsica_spanish_chiSquared3_df_0.99_Z, "corsica_spanish_chiSquared3_df_0.99_Z.csv")
  
  corsica_house_chiSquared3_df_0.99_autosomes <- chi_squared_test3_house_calc(corsica_house_alleleCounts_0.99_autosomes)
  write.csv(corsica_house_chiSquared3_df_0.99_autosomes, "corsica_house_chiSquared3_df_0.99_autosomes.csv")
  
  corsica_house_chiSquared3_df_0.99_Z <- chi_squared_test3_house_calc(corsica_house_alleleCounts_0.99_Z)
  write.csv(corsica_house_chiSquared3_df_0.99_Z, "corsica_house_chiSquared3_df_0.99_Z.csv")
  
  
  crete_spanish_chiSquared3_df_0.99_autosomes <- chi_squared_test3_spanish_calc(crete_spanish_alleleCounts_0.99_autosomes)
  write.csv(crete_spanish_chiSquared3_df_0.99_autosomes, "crete_spanish_chiSquared3_df_0.99_autosomes.csv")
  
  crete_spanish_chiSquared3_df_0.99_Z <- chi_squared_test3_spanish_calc(crete_spanish_alleleCounts_0.99_Z)
  write.csv(crete_spanish_chiSquared3_df_0.99_Z, "crete_spanish_chiSquared3_df_0.99_Z.csv")
  
  crete_house_chiSquared3_df_0.99_autosomes <- chi_squared_test3_house_calc(crete_house_alleleCounts_0.99_autosomes)
  write.csv(crete_house_chiSquared3_df_0.99_autosomes, "crete_house_chiSquared3_df_0.99_autosomes.csv")
  
  crete_house_chiSquared3_df_0.99_Z <- chi_squared_test3_house_calc(crete_house_alleleCounts_0.99_Z)
  write.csv(crete_house_chiSquared3_df_0.99_Z, "crete_house_chiSquared3_df_0.99_Z.csv")
  
  
  malta_spanish_chiSquared3_df_0.99_autosomes <- chi_squared_test3_spanish_calc(malta_spanish_alleleCounts_0.99_autosomes)
  write.csv(malta_spanish_chiSquared3_df_0.99_autosomes, "malta_spanish_chiSquared3_df_0.99_autosomes.csv")
  
  malta_spanish_chiSquared3_df_0.99_Z <- chi_squared_test3_spanish_calc(malta_spanish_alleleCounts_0.99_Z)
  write.csv(malta_spanish_chiSquared3_df_0.99_Z, "malta_spanish_chiSquared3_df_0.99_Z")
  
  malta_house_chiSquared3_df_0.99_autosomes <- chi_squared_test3_house_calc(malta_house_alleleCounts_0.99_autosomes)
  write.csv(malta_house_chiSquared3_df_0.99_autosomes, "malta_house_chiSquared3_df_0.99_autosomes.csv")
  
  malta_house_chiSquared3_df_0.99_Z <- chi_squared_test3_house_calc(malta_house_alleleCounts_0.99_Z)
  write.csv(malta_house_chiSquared3_df_0.99_Z, "malta_house_chiSquared3_df_0.99_Z.csv")
  
  
  sicily_spanish_chiSquared3_df_0.99_autosomes <- chi_squared_test3_spanish_calc(sicily_spanish_alleleCounts_0.99_autosomes)
  write.csv(sicily_spanish_chiSquared3_df_0.99_autosomes, "sicily_spanish_chiSquared3_df_0.99_autosomes.csv")
  
  sicily_spanish_chiSquared3_df_0.99_Z <- chi_squared_test3_spanish_calc(sicily_spanish_alleleCounts_0.99_Z)
  write.csv(sicily_spanish_chiSquared3_df_0.99_Z, "sicily_spanish_chiSquared3_df_0.99_Z.csv")
  
  sicily_house_chiSquared3_df_0.99_autosomes <- chi_squared_test3_house_calc(sicily_house_alleleCounts_0.99_autosomes)
  write.csv(sicily_house_chiSquared3_df_0.99_autosomes, "sicily_house_chiSquared3_df_0.99_autosomes.csv")
  
  sicily_house_chiSquared3_df_0.99_Z <- chi_squared_test3_house_calc(sicily_house_alleleCounts_0.99_Z)
  write.csv(sicily_house_chiSquared3_df_0.99_Z, "sicily_house_chiSquared3_df_0.99_Z.csv")
  
  
  # Create contigency dataframes, sep by parent - for built-in calc ####
  
  # From manual: "If x is a matrix with one row or column, or if x is a vector and y is not given, then a goodness-of-fit test is performed"
  
  # Function to create contignency "table" for Spanish parentage
  spanish_contigency_test3_df <- function(allele_count_df) {
    
    count_matrix <- matrix(c(allele_count_df$SA_total, allele_count_df$SD_total), nrow = 1, byrow = TRUE,
                           dimnames = list(c("spanish"),
                                           c("ancestral", "derived")))
    return(count_matrix)
  }
  
  # Function to create contignency "table" for house parentage
  house_contigency_test3_df <- function(allele_count_df) {
    
    count_matrix <- matrix(c(allele_count_df$HA_total, allele_count_df$HD_total), nrow = 1, byrow = TRUE,
                           dimnames = list(c("house"),
                                           c("ancestral", "derived")))
    return(count_matrix)
  }
  
  # Run the functions to create contingency tables
  corsica_spanish_contigency3_df_0.99 <- spanish_contigency_test3_df(corsica_alleleCounts_0.99)
  corsica_spanish_contigency3_df_0.99_autosomes <- spanish_contigency_test3_df(corsica_alleleCounts_0.99_autosomes)
  corsica_spanish_contigency3_df_0.99_Z <- spanish_contigency_test3_df(corsica_alleleCounts_0.99_Z)
  
  crete_spanish_contigency3_df_0.99 <- spanish_contigency_test3_df(crete_alleleCounts_0.99)
  crete_spanish_contigency3_df_0.99_autosomes <- spanish_contigency_test3_df(crete_alleleCounts_0.99_autosomes)
  crete_spanish_contigency3_df_0.99_Z <- spanish_contigency_test3_df(crete_alleleCounts_0.99_Z)
  
  malta_spanish_contigency3_df_0.99 <- spanish_contigency_test3_df(malta_alleleCounts_0.99)
  malta_spanish_contigency3_df_0.99_autosomes <- spanish_contigency_test3_df(malta_alleleCounts_0.99_autosomes)
  malta_spanish_contigency3_df_0.99_Z <- spanish_contigency_test3_df(malta_alleleCounts_0.99_Z)
  
  sicily_spanish_contigency3_df_0.99 <- spanish_contigency_test3_df(sicily_alleleCounts_0.99)
  sicily_spanish_contigency3_df_0.99_autosomes <- spanish_contigency_test3_df(sicily_alleleCounts_0.99_autosomes)
  sicily_spanish_contigency3_df_0.99_Z <- spanish_contigency_test3_df(sicily_alleleCounts_0.99_Z)
  
  corsica_house_contigency3_df_0.99 <- house_contigency_test3_df(corsica_alleleCounts_0.99)
  corsica_house_contigency3_df_0.99_autosomes <- house_contigency_test3_df(corsica_alleleCounts_0.99_autosomes)
  corsica_house_contigency3_df_0.99_Z <- house_contigency_test3_df(corsica_alleleCounts_0.99_Z)
  
  crete_house_contigency3_df_0.99 <- house_contigency_test3_df(crete_alleleCounts_0.99)
  crete_house_contigency3_df_0.99_autosomes <- house_contigency_test3_df(crete_alleleCounts_0.99_autosomes)
  crete_house_contigency3_df_0.99_Z <- house_contigency_test3_df(crete_alleleCounts_0.99_Z)
  
  malta_house_contigency3_df_0.99 <- house_contigency_test3_df(malta_alleleCounts_0.99)
  malta_house_contigency3_df_0.99_autosomes <- house_contigency_test3_df(malta_alleleCounts_0.99_autosomes)
  malta_house_contigency3_df_0.99_Z <- house_contigency_test3_df(malta_alleleCounts_0.99_Z)
  
  sicily_house_contigency3_df_0.99 <- house_contigency_test3_df(sicily_alleleCounts_0.99)
  sicily_house_contigency3_df_0.99_autosomes <- house_contigency_test3_df(sicily_alleleCounts_0.99_autosomes)
  sicily_house_contigency3_df_0.99_Z <- house_contigency_test3_df(sicily_alleleCounts_0.99_Z)
  
  # Run chisq.test & save output - for built-in calc ####
  # This function uses contingency tables to calc chi-sq
  
  sink("corsica_spanish_chisq3_auto_0.99.out")
  corsica_spanish_chisq_autosomes <- chisq.test(corsica_spanish_contigency3_df_0.99_autosomes)
  print(corsica_spanish_chisq_autosomes)
  sink()
  
  sink("corsica_spanish_chisq3_Z_0.99.out")
  corsica_spanish_chisq_Z <- chisq.test(corsica_spanish_contigency3_df_0.99_Z)
  print(corsica_spanish_chisq_Z)
  sink()
  
  sink("crete_spanish_chisq3_auto_0.99.out")
  crete_spanish_chisq_autosomes <- chisq.test(crete_spanish_contigency3_df_0.99_autosomes)
  print(crete_spanish_chisq_autosomes)
  sink()
  
  sink("crete_spanish_chisq3_Z_0.99.out")
  crete_spanish_chisq_Z <- chisq.test(crete_spanish_contigency3_df_0.99_Z)
  print(crete_spanish_chisq_Z)
  sink()
  
  sink("malta_spanish_chisq3_auto_0.99.out")
  malta_spanish_chisq_autosomes <- chisq.test(malta_spanish_contigency3_df_0.99_autosomes)
  print(malta_spanish_chisq_autosomes)
  sink()
  
  sink("malta_spanish_chisq3_Z_0.99.out")
  malta_spanish_chisq_Z <- chisq.test(malta_spanish_contigency3_df_0.99_Z)
  print(malta_spanish_chisq_Z)
  sink()
  
  sink("sicily_spanish_chisq3_auto_0.99.out")
  sicily_spanish_chisq_autosomes <- chisq.test(sicily_spanish_contigency3_df_0.99_autosomes)
  print(sicily_spanish_chisq_autosomes)
  sink()
  
  sink("sicily_spanish_chisq3_Z_0.99.out")
  sicily_spanish_chisq_Z <- chisq.test(sicily_spanish_contigency3_df_0.99_Z)
  print(sicily_spanish_chisq_Z)
  sink()

  
  sink("corsica_house_chisq3_auto_0.99.out")
  corsica_house_chisq_autosomes <- chisq.test(corsica_house_contigency3_df_0.99_autosomes)
  print(corsica_house_chisq_autosomes)
  sink()
  
  sink("corsica_house_chisq3_Z_0.99.out")
  corsica_house_chisq_Z <- chisq.test(corsica_house_contigency3_df_0.99_Z)
  print(corsica_house_chisq_Z)
  sink()
  
  sink("crete_house_chisq3_auto_0.99.out")
  crete_house_chisq_autosomes <- chisq.test(crete_house_contigency3_df_0.99_autosomes)
  print(crete_house_chisq_autosomes)
  sink()
  
  sink("crete_house_chisq3_Z_0.99.out")
  crete_house_chisq_Z <- chisq.test(crete_house_contigency3_df_0.99_Z)
  print(crete_house_chisq_Z)
  sink()
  
  sink("malta_house_chisq3_auto_0.99.out")
  malta_house_chisq_autosomes <- chisq.test(malta_house_contigency3_df_0.99_autosomes)
  print(malta_house_chisq_autosomes)
  sink()
  
  sink("malta_house_chisq3_Z_0.99.out")
  malta_house_chisq_Z <- chisq.test(malta_house_contigency3_df_0.99_Z)
  print(malta_house_chisq_Z)
  sink()
  
  sink("sicily_house_chisq3_auto_0.99.out")
  sicily_house_chisq_autosomes <- chisq.test(sicily_house_contigency3_df_0.99_autosomes)
  print(sicily_house_chisq_autosomes)
  sink()
  
  sink("sicily_house_chisq3_Z_0.99.out")
  sicily_house_chisq_Z <- chisq.test(sicily_house_contigency3_df_0.99_Z)
  print(sicily_house_chisq_Z)
  sink()

# Create chi-sq plot for test #3! ####
  
  # Check if user specifed lp data or not, adjust plot limits accordingly ####
  if (toupper(opt$lpData) == "YES") {
    plotLimits = c(0, 3200)
    
  } else if (toupper(opt$lpData) == "NO") {
    plotLimits = c(0, 1794000)
    
   
  } else {
    stop("Looks like you didn't specify if the data you provided is linkage pruned! Please specify yes or no")
  }
  
  # Function to create the plot ####

  plot_chi_sepParent_obs_exp <- function(spanish_crete_chi_df_auto, spanish_crete_chi_df_Z,
                                         spanish_corsica_chi_df_auto, spanish_corsica_chi_df_Z,
                                         spanish_sicily_chi_df_auto, spanish_sicily_chi_df_Z,
                                         spanish_malta_chi_df_auto, spanish_malta_chi_df_Z,
                                         house_crete_chi_df_auto, house_crete_chi_df_Z,
                                         house_corsica_chi_df_auto, house_corsica_chi_df_Z,
                                         house_sicily_chi_df_auto, house_sicily_chi_df_Z,
                                         house_malta_chi_df_auto, house_malta_chi_df_Z) {
    
    # Nested function to transform chi-sq df to long format
    toLong <- function(chi_df, populationName) {
      chi_df_long <- chi_df %>% 
        select(-expected_prop, -chiSquared) %>% # remove the chi-sq columns
        slice(1:(n() - 2)) %>% # remove the last two rows, NA row and ancestral row
        pivot_longer(cols = c(expected_number, observed_number), names_to = "Allele_type", values_to = "Allele_count") %>% # long format
        mutate(population = populationName) # add a label to preserve the population when we combine
      return(chi_df_long)
    }
    
    spanish_crete_chi_df_auto_long <- toLong(spanish_crete_chi_df_auto, "Crete")
    spanish_crete_chi_df_Z_long <- toLong(spanish_crete_chi_df_Z, "Crete")
    house_crete_chi_df_auto_long <- toLong(house_crete_chi_df_auto, "Crete")
    house_crete_chi_df_Z_long <- toLong(house_crete_chi_df_Z, "Crete")
    
    spanish_corsica_chi_df_auto_long <- toLong(spanish_corsica_chi_df_auto, "Corsica")
    spanish_corsica_chi_df_Z_long <- toLong(spanish_corsica_chi_df_Z, "Corsica")
    house_corsica_chi_df_auto_long <- toLong(house_corsica_chi_df_auto, "Corsica")
    house_corsica_chi_df_Z_long <- toLong(house_corsica_chi_df_Z, "Corsica")
    
    spanish_sicily_chi_df_auto_long <- toLong(spanish_sicily_chi_df_auto, "Sicily")
    spanish_sicily_chi_df_Z_long <- toLong(spanish_sicily_chi_df_Z, "Sicily")
    house_sicily_chi_df_auto_long <- toLong(house_sicily_chi_df_auto, "Sicily")
    house_sicily_chi_df_Z_long <- toLong(house_sicily_chi_df_Z, "Sicily")
    
    spanish_malta_chi_df_auto_long <- toLong(spanish_malta_chi_df_auto, "Malta")
    spanish_malta_chi_df_Z_long <- toLong(spanish_malta_chi_df_Z, "Malta")
    house_malta_chi_df_auto_long <- toLong(house_malta_chi_df_auto, "Malta")
    house_malta_chi_df_Z_long <- toLong(house_malta_chi_df_Z, "Malta")

    
    # Combine the long form dataframes and then separate by HD and SD!
    combined_spanish_df_auto <- bind_rows(spanish_crete_chi_df_auto_long, spanish_corsica_chi_df_auto_long, spanish_sicily_chi_df_auto_long, spanish_malta_chi_df_auto_long) %>%
      mutate(population = factor(population, levels = c("Crete", "Corsica", "Sicily", "Malta")))
    
    combined_spanish_df_Z <- bind_rows(spanish_crete_chi_df_Z_long, spanish_corsica_chi_df_Z_long, spanish_sicily_chi_df_Z_long, spanish_malta_chi_df_Z_long) %>%
      mutate(population = factor(population, levels = c("Crete", "Corsica", "Sicily", "Malta")))
    
    combined_house_df_auto <- bind_rows(house_crete_chi_df_auto_long, house_corsica_chi_df_auto_long, house_sicily_chi_df_auto_long, house_malta_chi_df_auto_long) %>%
      mutate(population = factor(population, levels = c("Crete", "Corsica", "Sicily", "Malta")))
    
    combined_house_df_Z <- bind_rows(house_crete_chi_df_Z_long, house_corsica_chi_df_Z_long, house_sicily_chi_df_Z_long, house_malta_chi_df_Z_long) %>%
      mutate(population = factor(population, levels = c("Crete", "Corsica", "Sicily", "Malta")))
    
    
    # Set color scheme for observed and expected for each population
    population_colors_expected <- c("Corsica" = "#87B8D1", "Crete" = "#6587A8", "Malta" = "#B46330", "Sicily" = "#D8B56B")
    population_colors_observed <- c("Corsica" = "#BBD8E6", "Crete" = "#9FBFD5", "Malta" = "#ED9668", "Sicily" = "#FAE2A0")
    
    # Create 4 separate plots & merge
    autoPlot_HD <- ggplot(combined_house_df_auto, aes(x = combination, y = Allele_count, color = population, fill = population, group = population)) +
      geom_bar(data = subset(combined_house_df_auto, Allele_type == "observed_number"), stat = "identity", position = "dodge", show.legend = FALSE, width = 0.9) + # no legend
      geom_point(data = subset(combined_house_df_auto, Allele_type == "expected_number"), aes(color = population), size = 3, position = position_dodge(width = 0.9), show.legend = FALSE) +
      labs(x = "", y = "") +
      scale_fill_manual(values = population_colors_observed) +  # Bar (observed) colors
      scale_color_manual(values = population_colors_expected) + # Point (expected) colots
      scale_x_discrete(labels = c("House_derived" = "Autosomes")) +
      scale_y_continuous(labels = scales::comma_format(), limits = plotLimits) +  # Setting y-axis limits
      #scale_y_continuous(labels = scales::comma_format()) +
      theme_minimal() +
      theme(axis.text.x = element_text(size = 14, color = "#4B68D9", face = "bold")) # increase axis label font size
    
    zPlot_HD <- ggplot(combined_house_df_Z, aes(x = combination, y = Allele_count, color = population, fill = population, group = population)) +
      geom_bar(data = subset(combined_house_df_Z, Allele_type == "observed_number"), stat = "identity", position = "dodge", show.legend = FALSE, width = 0.9) + # no legend
      geom_point(data = subset(combined_house_df_Z, Allele_type == "expected_number"), aes(color = population), size = 3, position = position_dodge(width = 0.9), show.legend = FALSE) +
      labs(x = "", y = "") + # no y label for this one!
      scale_fill_manual(values = population_colors_observed) +  # Bar (observed) colors
      scale_color_manual(values = population_colors_expected) + # Point (expected) colots
      scale_x_discrete(labels = c("House_derived" = "Z")) +
      scale_y_continuous(labels = scales::comma_format(), limits = plotLimits) +  # Remove numerical labels
      #scale_y_continuous(labels = scales::comma_format()) +
      theme_minimal() +
      theme(axis.text.x = element_text(size = 14, color = "#4B68D9", face = "bold"),
            axis.text.y = element_blank())
    
    autoPlot_SD <- ggplot(combined_spanish_df_auto, aes(x = combination, y = Allele_count, color = population, fill = population, group = population)) +
      geom_bar(data = subset(combined_spanish_df_auto, Allele_type == "observed_number"), stat = "identity", position = "dodge", show.legend = FALSE, width = 0.9,) + # no legend
      geom_point(data = subset(combined_spanish_df_auto, Allele_type == "expected_number"), aes(color = population), size = 3, position = position_dodge(width = 0.9), show.legend = FALSE) +
      labs(x = "", y = "") +
      scale_fill_manual(values = population_colors_observed) +  # Bar (observed) colors
      scale_color_manual(values = population_colors_expected) + # Point (expected) colots
      scale_x_discrete(labels = c("Spanish_derived" = "Autosomes")) +
      scale_y_continuous(labels = scales::comma_format(), limits = plotLimits) +  # Remove numerical labels
      #scale_y_continuous(labels = scales::comma_format()) +
      theme_minimal() +
      theme(axis.text.x = element_text(size = 14, color = "#BE441E", face = "bold"),
            axis.text.y = element_blank())
    
    zPlot_SD <- ggplot(combined_spanish_df_Z, aes(x = combination, y = Allele_count, color = population, fill = population, group = population)) +
      geom_bar(data = subset(combined_spanish_df_Z, Allele_type == "observed_number"), stat = "identity", position = "dodge", width = 0.9) + # keep legend, instruct ggarange to use shared legend to not affect bar width!
      geom_point(data = subset(combined_spanish_df_Z, Allele_type == "expected_number"), aes(color = population), size = 3, position = position_dodge(width = 0.9)) +
      labs(x = "", y = "", color = "", fill = "") + # no axis labels or legend title
      scale_fill_manual(values = population_colors_observed) +  # Bar (observed) colors
      scale_color_manual(values = population_colors_expected) + # Point (expected) colots
      scale_x_discrete(labels = c("Spanish_derived" = "Z")) +
      scale_y_continuous(labels = scales::comma_format(), limits = plotLimits) +  # Remove numerical labels
      #scale_y_continuous(labels = scales::comma_format()) +
      theme_minimal() +
      theme(axis.text.x = element_text(size = 14, color = "#BE441E", face = "bold"),
            axis.text.y = element_blank(),
            legend.text = element_text(size = 12))
    
    
    comboPlot <- ggarrange(autoPlot_HD, NULL, zPlot_HD, NULL, autoPlot_SD, NULL, zPlot_SD,
                           widths = c(1.15,-0.2, 1, -0.2, 1, -0.2, 1),
                           common.legend = TRUE, legend = "right",
                           ncol = 7, nrow = 1)
    
    comboPlot <- annotate_figure(comboPlot,
                                 left = text_grob("Allele Count", color = "black", rot = 90, size = 14)
    )
    
    return(comboPlot)
    
  }
  
  # Run the function & save output ####
  png(filename = "chisq3_plot.png", width = 3500, height = 2000, res = 300)
  plot_chi_sepParent_obs_exp(crete_spanish_chiSquared3_df_0.99_autosomes, crete_spanish_chiSquared3_df_0.99_Z,
                             corsica_spanish_chiSquared3_df_0.99_autosomes, corsica_spanish_chiSquared3_df_0.99_Z,
                             sicily_spanish_chiSquared3_df_0.99_autosomes, sicily_spanish_chiSquared3_df_0.99_Z,
                             malta_spanish_chiSquared3_df_0.99_autosomes, malta_spanish_chiSquared3_df_0.99_Z,
                             crete_house_chiSquared3_df_0.99_autosomes, crete_house_chiSquared3_df_0.99_Z, 
                             corsica_house_chiSquared3_df_0.99_autosomes, corsica_house_chiSquared3_df_0.99_Z,
                             sicily_house_chiSquared3_df_0.99_autosomes, sicily_house_chiSquared3_df_0.99_Z,
                             malta_house_chiSquared3_df_0.99_autosomes, malta_house_chiSquared3_df_0.99_Z)
  dev.off()
  