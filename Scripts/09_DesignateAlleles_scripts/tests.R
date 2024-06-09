# Load packages ####
library(tidyverse) # data frame manipulation & plotting, version 2.0.0
library(stats) # one sample prop test functions, version 4.3.3
library(gplots) # balloonplot, version 3.1.3.1
library(graphics) # mosaicplot, version 4.3.3
library(rcompanion) #cramersv, version 2.4.36
library(optparse) # allow use of command line flags

# Run from command line! ####

flag_list <- list(
  
  make_option(c("--sansSAHA_0.90"), type = "character", default = NULL, 
              help = "Input designation file for 0.90 cutoff, sansSAHA", metavar = "file_path"),
  
  make_option(c("--SAHA_0.90"), type = "character", default = NULL, 
              help = "Input designation file for 0.90 cutoff, SAHA", metavar = "file_path"),
  
  make_option(c("--sansSAHA_0.99"), type = "character", default = NULL, 
              help = "Input designation file for 0.99 cutoff, sansSAHA", metavar = "file_path"),
  
  make_option(c("--SAHA_0.99"), type = "character", default = NULL, 
              help = "Input designation file for 0.99 cutoff, SAHA", metavar = "file_path"),
  
  make_option(c("--width"), type = "integer", default = 1600, 
              help = "Width of output images in pixels [default= %default]", metavar = "integer"),
  
  make_option(c("--height"), type = "integer", default = 1350, 
              help = "Height of output images in pixels [default= %default]", metavar = "integer"),
  
  make_option(c("--res"), type = "integer", default = 300, 
              help = "Resolution of output images in DPI [default= %default]", metavar = "integer")
)

# Set the flag list
opt_parser <- OptionParser(option_list = flag_list)
opt <- parse_args(opt_parser) # can reference user provided arguments using opt$res for example!


# Import Data ####

# Create a function for reading in designation files
  read_desFile <- function(file_path) {
    dat <- read.csv(file_path, header = TRUE, sep = "\t") # will replace the hash before Locus with X.
    dat <- dat %>%
      rename(
        "locus" = "X.Locus", # rename column to just "locus"
        "K034" = "X034" # rename column to "K034" (the K was dropped by vcftools when generating VCFs)
      )
    return(dat)
  }

# Read in the designation files
  cutoff_0.90_sansSAHA_dat <- read_desFile(opt$sansSAHA_0.90)
  
  cutoff_0.90_SAHA_dat <- read_desFile(opt$SAHA_0.90)
  
  cutoff_0.99_sansSAHA_dat <- read_desFile(opt$sansSAHA_0.99)
  
  cutoff_0.99_SAHA_dat <- read_desFile(opt$SAHA_0.99)

# Chi Squared ####
  
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
  
  # Run the function for each population at each cutoff
  corsica_alleleCounts_0.99 <- count_Italian_alleles(cutoff_0.99_SAHA_dat, population_prefix = "K0")
  corsica_alleleCounts_0.90 <- count_Italian_alleles(cutoff_0.90_SAHA_dat, population_prefix = "K0")
  
  crete_alleleCounts_0.99 <- count_Italian_alleles(cutoff_0.99_SAHA_dat, population_prefix = "C0")
  crete_alleleCounts_0.90 <- count_Italian_alleles(cutoff_0.90_SAHA_dat, population_prefix = "C0")
  
  malta_alleleCounts_0.99 <- count_Italian_alleles(cutoff_0.99_SAHA_dat, population_prefix = "M0")
  malta_alleleCounts_0.90 <- count_Italian_alleles(cutoff_0.90_SAHA_dat, population_prefix = "M0")
  
  sicily_alleleCounts_0.99 <- count_Italian_alleles(cutoff_0.99_SAHA_dat, population_prefix = "S0")
  sicily_alleleCounts_0.90 <- count_Italian_alleles(cutoff_0.90_SAHA_dat, population_prefix = "S0")
  
  # First, compute chi-squared manually! ####
  # Function to create empty chi_squared df dataframes for each population/cutoff combo
  
  establish_chiSquare <- function() {
    population_chiSquared_df_cutoff <- data.frame(
      combination = c("House_derived", "House_ancestral", "Spanish_derived", "Spanish_ancestral")
    )
    return(population_chiSquared_df_cutoff)
  }

  corsica_chiSquared_df_0.99 <- establish_chiSquare()
  corsica_chiSquared_df_0.90 <- establish_chiSquare()
  
  crete_chiSquared_df_0.99 <- establish_chiSquare()
  crete_chiSquared_df_0.90 <- establish_chiSquare()
  
  malta_chiSquared_df_0.99 <- establish_chiSquare()
  malta_chiSquared_df_0.90 <- establish_chiSquare()
  
  sicily_chiSquared_df_0.99 <- establish_chiSquare()
  sicily_chiSquared_df_0.90 <- establish_chiSquare()
  
  # Function to create chi squared table!
  propDerived <- function(chi_squared_df, alleleCounts_df) {
    
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
  
  corsica_chiSquared_df_0.99 <- propDerived(chi_squared_df = corsica_chiSquared_df_0.99, alleleCounts_df = corsica_alleleCounts_0.99)
  corsica_chiSquared_df_0.90 <- propDerived(chi_squared_df = corsica_chiSquared_df_0.90, alleleCounts_df = corsica_alleleCounts_0.90)
  
  crete_chiSquared_df_0.99 <- propDerived(chi_squared_df = crete_chiSquared_df_0.99, alleleCounts_df = crete_alleleCounts_0.99)
  crete_chiSquared_df_0.90 <- propDerived(chi_squared_df = crete_chiSquared_df_0.90, alleleCounts_df = crete_alleleCounts_0.90)
  
  malta_chiSquared_df_0.99 <- propDerived(chi_squared_df = malta_chiSquared_df_0.99, alleleCounts_df = malta_alleleCounts_0.99)
  malta_chiSquared_df_0.90 <- propDerived(chi_squared_df = malta_chiSquared_df_0.90, alleleCounts_df = malta_alleleCounts_0.90)
  
  sicily_chiSquared_df_0.99 <- propDerived(chi_squared_df = sicily_chiSquared_df_0.99, alleleCounts_df = sicily_alleleCounts_0.99)
  sicily_chiSquared_df_0.90 <- propDerived(chi_squared_df = sicily_chiSquared_df_0.90, alleleCounts_df = sicily_alleleCounts_0.90)
  
  # Write the output to tsv files!
  write.table(corsica_chiSquared_df_0.99, file="corsica_manual_chiSquared_0.99.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(corsica_chiSquared_df_0.90, file="corsica_manual_chiSquared_0.90.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  write.table(crete_chiSquared_df_0.99, file="crete_manual_chiSquared_0.99.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(crete_chiSquared_df_0.90, file="crete_manual_chiSquared_0.90.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  write.table(malta_chiSquared_df_0.99, file="malta_manual_chiSquared_0.99.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(malta_chiSquared_df_0.90, file="malta_manual_chiSquared_0.90.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  write.table(sicily_chiSquared_df_0.99, file="sicily_manual_chiSquared_0.99.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(sicily_chiSquared_df_0.90, file="sicily_manual_chiSquared_0.90.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  
  # Compute chi-squared using built in packages ####
  
  # Function to create contignency "tables" (actually using df format)
  # Reference: http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
  contigency_df <-function(population_alleleCounts_cutoff) {
    
    count_matrix <- matrix(c(population_alleleCounts_cutoff$HA_total, population_alleleCounts_cutoff$HD_total, population_alleleCounts_cutoff$SA_total, population_alleleCounts_cutoff$SD_total), nrow = 2, byrow = TRUE,
                             dimnames = list(c("house", "spanish"),
                                             c("ancestral", "derived")))
    return(count_matrix)
  }
  
  # Run the function
  corsica_contigency_0.99 <- contigency_df(corsica_alleleCounts_0.99)
  corsica_contigency_0.90 <- contigency_df(corsica_alleleCounts_0.90)
  
  crete_contigency_0.99 <- contigency_df(crete_alleleCounts_0.99)
  crete_contigency_0.90 <- contigency_df(crete_alleleCounts_0.90)
  
  malta_contigency_0.99 <- contigency_df(malta_alleleCounts_0.99)
  malta_contigency_0.90 <- contigency_df(malta_alleleCounts_0.90)
  
  sicily_contigency_0.99 <- contigency_df(sicily_alleleCounts_0.99)
  sicily_contigency_0.90 <- contigency_df(sicily_alleleCounts_0.90)
  
  # Visualize the contigency tables!
  
  visualize <- function(contigency_matrix, plot_title, baloon_fileName, mosaic_fileName) { # titles should be string in ""
    dt <- as.table(contigency_matrix)
    
    # Plot & save!
    png(filename = baloon_fileName, width = opt$width, height = opt$height, res = opt$res)
    balloonplot(t(dt), main = plot_title, xlab ="", ylab="",
                label = FALSE, show.margins = FALSE)
    dev.off()
    
    png(filename = mosaic_fileName, width = opt$width, height = opt$height, res = opt$res)
    mosaicplot(dt, shade = TRUE, las=2,
               main = plot_title)
    dev.off()
    
  }

  visualize(corsica_contigency_0.99, plot_title = "Corsica (0.99 cutoff) - allele designations",
            baloon_fileName = "corsica_0.99_chiContigencyTable_baloon.png",
            mosaic_fileName = "corsica_0.99_chiContigencyTable_mosaic.png")
  
  visualize(corsica_contigency_0.90, plot_title = "Corsica (0.90 cutoff) - allele designations",
            baloon_fileName = "corsica_0.90_chiContigencyTable_baloon.png",
            mosaic_fileName = "corsica_0.90_chiContigencyTable_mosaic.png")
  
  visualize(crete_contigency_0.99, plot_title = "Crete (0.99 cutoff) - allele designations",
            baloon_fileName = "crete_0.99_chiContigencyTable_baloon.png",
            mosaic_fileName = "crete_0.99_chiContigencyTable_mosaic.png")
  
  visualize(crete_contigency_0.90, plot_title = "Crete (0.90 cutoff) - allele designations",
            baloon_fileName = "crete_0.90_chiContigencyTable_baloon.png",
            mosaic_fileName = "crete_0.90_chiContigencyTable_mosaic.png")
  
  visualize(malta_contigency_0.99, plot_title = "Malta (0.99 cutoff) - allele designations",
            baloon_fileName = "malta_0.99_chiContigencyTable_baloon.png",
            mosaic_fileName = "malta_0.99_chiContigencyTable_mosaic.png")
  
  visualize(malta_contigency_0.90, plot_title = "Malta (0.90 cutoff) - allele designations",
            baloon_fileName = "malta_0.90_chiContigencyTable_baloon.png",
            mosaic_fileName = "malta_0.90_chiContigencyTable_mosaic.png")
  
  visualize(sicily_contigency_0.99, plot_title = "Sicily (0.99 cutoff) - allele designations",
            baloon_fileName = "sicily_0.99_chiContigencyTable_baloon.png",
            mosaic_fileName = "sicily_0.99_chiContigencyTable_mosaic.png")
  
  visualize(sicily_contigency_0.90, plot_title = "Sicily (0.90 cutoff) - allele designations",
            baloon_fileName = "sicily_0.90_chiContigencyTable_baloon.png",
            mosaic_fileName = "sicily_0.90_chiContigencyTable_mosaic.png")
  
  # Run chi-squared tests, Pearson's Chi-squared test with Yates' continuity correction & save output!
  
  sink("corsica_0.99_chisq.test.out")
  corsica_0.99_chisq <- chisq.test(corsica_contigency_0.99) # X-squared = 76534, df = 1, p-value < 2.2e-16
  print(corsica_0.99_chisq)
  sink()
  
  sink("corsica_0.90_chisq.test.out")
  corsica_0.90_chisq <- chisq.test(corsica_contigency_0.90) # X-squared = 219825, df = 1, p-value < 2.2e-16
  print(corsica_0.90_chisq)
  sink()
  
  sink("crete_0.99_chisq.test.out")
  crete_0.99_chisq <- chisq.test(crete_contigency_0.99) # X-squared = 15610, df = 1, p-value < 2.2e-16
  print(crete_0.99_chisq)
  sink()
  
  sink("crete_0.90_chisq.test.out")
  crete_0.90_chisq <- chisq.test(crete_contigency_0.90) # X-squared = 68868, df = 1, p-value < 2.2e-16
  print(crete_0.90_chisq)
  sink()
  
  sink("malta_0.99_chisq.test.out")
  malta_0.99_chisq <- chisq.test(malta_contigency_0.99) # X-squared = 43406, df = 1, p-value < 2.2e-16
  print(malta_0.99_chisq)
  sink()
  
  sink("malta_0.90_chisq.test.out")
  malta_0.90_chisq <- chisq.test(malta_contigency_0.90) # X-squared = 116713, df = 1, p-value < 2.2e-16
  print(malta_0.90_chisq)
  sink()

  sink("sicily_0.99_chisq.test.out")
  sicily_0.99_chisq <- chisq.test(sicily_contigency_0.99) # X-squared = 40556, df = 1, p-value < 2.2e-16
  print(sicily_0.99_chisq)
  sink()
  
  sink("sicily_0.90_chisq.test.out")
  sicily_0.90_chisq <- chisq.test(sicily_contigency_0.90) # X-squared = 112549, df = 1, p-value < 2.2e-16
  print(sicily_0.90_chisq)
  sink()
  
# One sample proportion test (test derived against expectation of 50%) ####
  
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
  corsica_DA_Counts_0.90 <- count_Italian_DA(cutoff_0.90_SAHA_dat, population_prefix = "K0")
  
  crete_DA_Counts_0.99 <- count_Italian_DA(cutoff_0.99_SAHA_dat, population_prefix = "C0")
  crete_DA_Counts_0.90 <- count_Italian_DA(cutoff_0.90_SAHA_dat, population_prefix = "C0")
  
  malta_DA_Counts_0.99 <- count_Italian_DA(cutoff_0.99_SAHA_dat, population_prefix = "M0")
  malta_DA_Counts_0.90 <- count_Italian_DA(cutoff_0.90_SAHA_dat, population_prefix = "M0")
  
  sicily_DA_Counts_0.99 <- count_Italian_DA(cutoff_0.99_SAHA_dat, population_prefix = "S0")
  sicily_DA_Counts_0.90 <- count_Italian_DA(cutoff_0.90_SAHA_dat, population_prefix = "S0")
  
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
  corsica_prop_0.99 <- prop_test(corsica_DA_Counts_0.99) # reject null (supports purging of derived), p ≈ 0.4769133, p-value: < 2.2e-16
  print(corsica_prop_0.99)
  sink()
  
  sink("corsica_0.90_prop.test.out")
  corsica_prop_0.90 <- prop_test(corsica_DA_Counts_0.90) # reject null (supports purging of derived), p ≈ 0.4754222, p-value: < 2.2e-16
  print(corsica_prop_0.90)
  sink()
  
  sink("crete_0.99_prop.test.out")
  crete_prop_0.99 <- prop_test(crete_DA_Counts_0.99) # reject null (supports purging of derived), p ≈ 0.45231, p-value: < 2.2e-16
  print(crete_prop_0.99)
  sink()
  
  sink("crete_0.90_prop.test.out")
  crete_prop_0.90 <- prop_test(crete_DA_Counts_0.90) # reject null (supports purging of derived), p ≈ 0.4549287, p-value: < 2.2e-16
  print(crete_prop_0.90)
  sink()
  
  sink("malta_0.99_prop.test.out")
  malta_prop_0.99 <- prop_test(malta_DA_Counts_0.99) # Fail to reject null, p ≈ 0.5069567, p-value: < 2.2e-16
  print(malta_prop_0.99)
  sink()
  
  sink("malta_0.90_prop.test.out")
  malta_prop_0.90 <- prop_test(malta_DA_Counts_0.90) # Fail to reject null, p ≈ 0.511139, p-value: < 2.2e-16
  print(malta_prop_0.90)
  sink()
  
  sink("sicily_0.99_prop.test.out")
  sicily_prop_0.99 <- prop_test(sicily_DA_Counts_0.99) # Fail to reject null, p ≈ 0.5150267, p-value: < 2.2e-16
  print(sicily_prop_0.99)
  sink()
  
  sink("sicily_0.90_prop.test.out")
  sicily_prop_0.90 <- prop_test(sicily_DA_Counts_0.90) # Fail to reject null, p ≈ 0.5172192, p-value: < 2.2e-16
  print(sicily_prop_0.90)
  sink()

# Hardy-Weinberg ####
  
  # Function to create chi-sq HW tables (possibility- DD, DA, AA)
  population_genotype_counts <- function(dat, population_prefix) {     # population_prefix should be entered as a string in ""

    observed_genotype_counts <- dat %>% # dat being the designation dataframe
      select(starts_with(population_prefix)) %>% # only pull columns from the designation df which start with the provided population prefix
      pivot_longer(cols = everything(), names_to = "sample", values_to = "genotype") %>%
      group_by(sample) %>% # group by sample, (e.g. corsica 0.99 - 2,885,631 rows at this point bc 21 samples * 137,411 loci)
      filter(genotype != "-") %>% # remove missing genotypes, (e.g. corsica 0.99 - 2,878,328 rows at this point)
      group_by(genotype) %>% # group by population and genotype
      summarize(count = n()) %>% # then count occurrence of each genotype for each population
      pivot_wider(names_from = genotype, values_from = count, values_fill = list(count = 0)) %>%
      mutate(DD = `HD/HD` + `SD/SD`,
             DA = `A/HD` + `A/SD` + `HD/A` + `SD/A`,
             AA = `A/A`) %>%
      mutate(genotype_total = `DD` +`DA` + `AA`) %>%
      select(DD, DA, AA, genotype_total)
    
    allele_freq_df <- observed_genotype_counts %>%
      mutate(D_count = ((DD * 2) + (DA * 1)), # create a total D allele count column
             A_count = ((AA * 2) + (DA * 1))) %>% # create a total A allele count column
      mutate(D_freq = (D_count/(genotype_total * 2)), # calculate allele freqs based on the counts
             A_freq = (A_count/(genotype_total * 2))) %>%
      select(D_freq, A_freq) # only keep the frequency columns, should add up to 1!
    
    # Create a one row dataframe with the expected genotype counts
    expected_genotype_counts_row <- data.frame(
      DD = ((allele_freq_df$D_freq)^2)*(observed_genotype_counts$genotype_total),
      DA = (2*allele_freq_df$D_freq*allele_freq_df$A_freq)*(observed_genotype_counts$genotype_total),
      AA = ((allele_freq_df$A_freq)^2)*(observed_genotype_counts$genotype_total)
    )
    
    # Append the expected genotype counts df to the observed genotype counts df
    OE_genotype_counts <- observed_genotype_counts %>%
      select(-genotype_total) %>% # remove genotype total column
      bind_rows(expected_genotype_counts_row) # adds row with observed counts
    
    # add rownames for easy calling in next step
    rownames(OE_genotype_counts) <- c("observed", "expected")
    
    # Calculate the chisq for each genotype in a one row dataframe
    chi_sq_row <- data.frame(
      DD = (((OE_genotype_counts["observed", "DD"] - OE_genotype_counts["expected", "DD"])^2)/OE_genotype_counts["expected", "DD"]),
      DA = (((OE_genotype_counts["observed", "DA"] - OE_genotype_counts["expected", "DA"])^2)/OE_genotype_counts["expected", "DA"]),
      AA = (((OE_genotype_counts["observed", "AA"] - OE_genotype_counts["expected", "AA"])^2)/OE_genotype_counts["expected", "AA"])
    ) %>%
      mutate(chisquared = DD + DA  + AA)
    
    chi_sq_df <- OE_genotype_counts %>%
      bind_rows(chi_sq_row)
    
    rownames(chi_sq_df) <- c("observed", "expected", "chisquared")
    
    return(chi_sq_df)
  } 
  
 
  # Get the chi-sq tables
  # df = (# rows -1) * (# cols -1)
  # df = (2-1) * (# 3-1) = 2
  corsica_HW_chisq_0.99 <- population_genotype_counts(cutoff_0.99_sansSAHA_dat, "K0") # 215,976
  corsica_HW_chisq_0.90 <- population_genotype_counts(cutoff_0.90_sansSAHA_dat, "K0") # 274,682

  crete_HW_chisq_0.99 <- population_genotype_counts(cutoff_0.99_sansSAHA_dat, "C0") # 593,366
  crete_HW_chisq_0.90 <- population_genotype_counts(cutoff_0.90_sansSAHA_dat, "C0") # 669,137
  
  malta_HW_chisq_0.99 <- population_genotype_counts(cutoff_0.99_sansSAHA_dat, "M0") # 364,741
  malta_HW_chisq_0.90 <- population_genotype_counts(cutoff_0.90_sansSAHA_dat, "M0") # 505,800
  
  sicily_HW_chisq_0.99 <- population_genotype_counts(cutoff_0.99_sansSAHA_dat, "S0") # 306,997
  sicily_HW_chisq_0.90 <- population_genotype_counts(cutoff_0.90_sansSAHA_dat, "S0") # 409,651
  
  # reject null, not in HWE (OBVIOUSLY)
  
  # Function that will create df showing observed and expected genotype frequencies
  genotype_freqs_observed_expected <- function(dat, population_prefix) {     # population_prefix should be entered as a string in ""
    
    observed_genotype_counts <- dat %>% # dat being the designation dataframe
      select(starts_with(population_prefix)) %>% # only pull columns from the designation df which start with the provided population prefix
      pivot_longer(cols = everything(), names_to = "sample", values_to = "genotype") %>%
      group_by(sample) %>% # group by sample, (e.g. corsica 0.99 - 2,885,631 rows at this point bc 21 samples * 137,411 loci)
      filter(genotype != "-") %>% # remove missing genotypes, (e.g. corsica 0.99 - 2,878,328 rows at this point)
      group_by(genotype) %>% # group by population and genotype
      summarize(count = n()) %>% # then count occurrence of each genotype for each population
      pivot_wider(names_from = genotype, values_from = count, values_fill = list(count = 0)) %>%
      mutate(DD = `HD/HD` + `SD/SD`,
             DA = `A/HD` + `A/SD` + `HD/A` + `SD/A`,
             AA = `A/A`) %>%
      mutate(genotype_total = `DD` +`DA` + `AA`) %>%
      select(DD, DA, AA, genotype_total)
    
    observed_genotype_freq <- observed_genotype_counts %>%
      mutate(DD_freq = (DD/genotype_total),
             DA_freq = (DA/genotype_total),
             AA_freq = (AA/genotype_total)) %>%
      select(DD_freq, DA_freq, AA_freq )
    
    allele_freq_df <- observed_genotype_counts %>%
      mutate(D_count = ((DD * 2) + (DA * 1)), # create a total D allele count column
             A_count = ((AA * 2) + (DA * 1))) %>% # create a total A allele count column
      mutate(D_freq = (D_count/(genotype_total * 2)), # calculate allele freqs based on the counts
             A_freq = (A_count/(genotype_total * 2))) %>%
      select(D_freq, A_freq) # only keep the frequency columns, should add up to 1!
    
    # Create a one row dataframe with the expected genotype counts
    expected_genotype_counts_row <- data.frame(
      DD = ((allele_freq_df$D_freq)^2)*(observed_genotype_counts$genotype_total), # expected count is allele freq squared (because 2 copies) times the genotype total
      DA = (2*allele_freq_df$D_freq*allele_freq_df$A_freq)*(observed_genotype_counts$genotype_total),
      AA = ((allele_freq_df$A_freq)^2)*(observed_genotype_counts$genotype_total)
    )
    
    expected_genotype_freq <- expected_genotype_counts_row %>%
      mutate(DD_freq = (DD/observed_genotype_counts$genotype_total), # expected freq is the expected count/total genotype
             DA_freq = (DA/observed_genotype_counts$genotype_total),
             AA_freq = (AA/observed_genotype_counts$genotype_total)) %>%
      select(DD_freq, DA_freq, AA_freq )
    
    observed_expected_freq <- observed_genotype_freq %>% # bind the observed and expected freq table together
      bind_rows(expected_genotype_freq)
    
    rownames(observed_expected_freq) <- c("observed", "expected") # add row names
    
    return(observed_expected_freq)
  }
  
  # Run the function
  corsica_genotypeFreq_0.99 <- genotype_freqs_observed_expected(cutoff_0.99_sansSAHA_dat, "K0") 
  corsica_genotypeFreq_0.90 <- genotype_freqs_observed_expected(cutoff_0.90_sansSAHA_dat, "K0") 
  
  crete_genotypeFreq_0.99 <- genotype_freqs_observed_expected(cutoff_0.99_sansSAHA_dat, "C0") 
  crete_genotypeFreq_0.90 <- genotype_freqs_observed_expected(cutoff_0.90_sansSAHA_dat, "C0")
  
  malta_genotypeFreq_0.99 <- genotype_freqs_observed_expected(cutoff_0.99_sansSAHA_dat, "M0")
  malta_genotypeFreq_0.90 <- genotype_freqs_observed_expected(cutoff_0.90_sansSAHA_dat, "M0")
  
  sicily_genotypeFreq_0.99 <- genotype_freqs_observed_expected(cutoff_0.99_sansSAHA_dat, "S0")
  sicily_genotypeFreq_0.90 <- genotype_freqs_observed_expected(cutoff_0.90_sansSAHA_dat, "S0")
  
  # Write the outputs to a tsv
  write.table(corsica_genotypeFreq_0.99, file="corsica_0.99_HW_genotypeFreq.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(corsica_genotypeFreq_0.90, file="corsica_0.90_HW_genotypeFreq.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  write.table(crete_genotypeFreq_0.99, file="crete_0.99_HW_genotypeFreq.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(crete_genotypeFreq_0.90, file="crete_0.90_HW_genotypeFreq.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  write.table(malta_genotypeFreq_0.99, file="malta_0.99_HW_genotypeFreq.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(malta_genotypeFreq_0.90, file="malta_0.90_HW_genotypeFreq.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  write.table(sicily_genotypeFreq_0.99, file="sicily_0.99_HW_genotypeFreq.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(sicily_genotypeFreq_0.90, file="sicily_0.90_HW_genotypeFreq.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
