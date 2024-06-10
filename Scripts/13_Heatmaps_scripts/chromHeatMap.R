#!/usr/bin/env Rscript

# R script that that is used to site designation across chromosome
# Use with pruned/unpruned designation .tsv files

# Example usage: Rscript chromHeatMap.R --input designationFiles/cutoff_0.90_SAHA.tsv 
# Install & load necessary packages ####

# List of packages needed
packages <- c("optparse", # v 1.7.5, allow use of command line flags
              "tidyverse" # v 2.0.0 tidyr, dplyr, purr (keep(~ nrow(.) > 1))
              ) 

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], repos = "https://cloud.r-project.org/") # error message if we don't set cran mirror
}

# Load the packages
invisible(lapply(packages, library, character.only = TRUE))

# Run from command line! ####

flag_list <- list(
  make_option(c("--input"), type = "character", default = NULL, 
              help = "Input designation file!", metavar = "file_path"),
  
  make_option(c("--width"), type = "integer", default = 6700, 
              help = "Width of output images in pixels [default= %default]", metavar = "integer"),
  
  make_option(c("--height"), type = "integer", default = 5000, 
              help = "Height of output images in pixels [default= %default]", metavar = "integer"),
  
  make_option(c("--res"), type = "integer", default = 300, 
              help = "Resolution of output images in DPI [default= %default]", metavar = "integer")
)

# Set the flag list
opt_parser <- OptionParser(option_list = flag_list)
opt <- parse_args(opt_parser)

# Imports ####

# Function to read in designation files!

read_desFile <- function(file_path) {
  dat <- read.csv(file_path, header = TRUE, sep = "\t") # will replace the hash before Locus with X.
  dat <- dat %>%
    rename(
      "locus" = "X.Locus", # rename column to just "locus"
      "K034" = "X034" # rename column to "K034" (the K was dropped by vcftools when generating VCFs)
    )
  
  return(dat)
}

# Read input
desFile_dat <- read_desFile(opt$input)

# Dataframe manipulation ####

# Function to pull chromosome, position, & site designation
chrom_pos_des <- function(dat){
  
  dat <- dat %>%
    select(locus, Site_Designation) # pull just locus and site designation
  
  dat <- dat %>% # separate the locus into its parts
    separate(locus, into = c("chromosome", "position", "ref", "alt"), sep = "_")
  
  dat <- dat %>% # we don't need the ref and alt allele so drop that
    select(-ref, -alt)
  
  dat$position <- as.numeric(dat$position)# we don't want the positions treated as characters
  
  return(dat)
}

# Run the function!
desFile_chrom_pos_des <- chrom_pos_des(desFile_dat)


# Split the dataframes by chromosomes/scaffold
print("Removing chromosomes that do not have more than 1 designated site")
split_by_chrom <- split(desFile_chrom_pos_des, desFile_chrom_pos_des$chromosome) %>% # creates a list which stores the dataframes
  keep(~ nrow(.) > 1) # Only keep dataframes with more than 1 row, meaning more than one posiiton has been designated


# Heatmap ####

# Create function
chromosome_heatmap <- function(chromosome_df) {
  
  count_df <- chromosome_df %>% 
    group_by(chromosome, position, Site_Designation) %>%
    summarise(count = n()) %>% # Will pull a count for the number of a designation type at each position (always 1)
    ungroup() # Only needed it grouped to count
  
  # Wide format! Chromosome, position, and count of HD and SD at that position. If HD count 1, then SD count 0 (obv)
  heatmap_data <- pivot_wider(count_df, names_from = Site_Designation, values_from = count, values_fill = 0)
  
  # Count matrix for positions and counts! Transpose so site designation is row
  heatmap_matrix <- t(as.matrix(heatmap_data[, -c(1, 2)])) # columns named V1, V2, etc. will correct next step
  
  # Pull the chrom position from chromosome_df and replace V1, V2, V3 as column names
  col_names <- as.character(chromosome_df$position) # character that stores positions e.g. c("614190", "17287392")
  colnames(heatmap_matrix) <- col_names # renames the columns of matrix
  
  # Check if heatmap matrix has at least two rows
  if (nrow(heatmap_matrix) < 2) { # aka need to have at least one SD site and one HD site
    warning("Skipping chromosome ", unique(chromosome_df$chromosome), "! Require at least one SD site and one HD site")
    return(NULL)
  } # without this, will halt function!
  
  # Rename the rows in the matrix! DON'T do this, SD and HD order changes in matrix based on which appears first in data
  #row_names <- c("HD/SA", "SD/HA") # if one allele (ref or alt) is House derived, then the other is Spanish ancestral and vice versa
  #rownames(heatmap_matrix) <- row_names
  
  # Set the margins, the chromosomes with less sites (aka positions) will be wayyyyy to zoomed in if not
  num_positions <- ncol(heatmap_matrix) # the number of columnns = number of positions
  margins <- if (num_positions > 16) c(6, 4) else c(12, 16) # trial and error, if 16 + positions set margins to c(6,4) if less set them to c(12,16)
  
  # Save the plot within the function!
  chromosome_name <- unique(chromosome_df$chromosome) # pull the chromosome name here so it doesn't cause issues in paste
  png_filename <- paste0("siteDes_heatmap_", chromosome_name, ".png")
  
  # Save the heatmap as a PNG
  png(filename = png_filename, width = opt$width, height = opt$height, res = opt$res)
  
  # Plot the heatmap within the PNG context
  heat_colors <- c("white", "blue")  
  heatmap(heatmap_matrix, col = heat_colors, Rowv = NA, Colv = NA,
          scale = "column", # positions are in columns!
          margins = margins,  
          main = paste("Site designations along", unique(chromosome_df$chromosome)), # without unique, will print the chromsome name per number of sites
          xlab = "Position on chromosome", ylab = "Site Designation",
          symm = FALSE # not a symmetric matrix
  )
  
  dev.off()  # Close the PNG device

}

heatmap_list <- map(split_by_chrom, chromosome_heatmap)


