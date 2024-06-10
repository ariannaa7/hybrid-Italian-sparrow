#!/usr/bin/env Rscript

# R script that uses Genotype Plot (https://github.com/JimWhiting91/genotype_plot?tab=readme-ov-file)
# To visualize genotypes along VCF. Can use .vcf or .vcf.gz inputs

# Example usage: 
# Rscript genotypePlot_VCFs.R --VCF VCFs_subsettedDerviedSites/cutoff_0.99_subsettedVCFs/CM071426.1_subset_cutoff_0.99.vcf


# To run on multiple VCFs in a directory: 
# for vcf in VCFs_subsettedDerviedSites/cutoff_0.99_subsettedVCFs/*vcf; do      Rscript genotypePlot_VCFs.R --VCF ${vcf};  done

# ERROR NOTE!!!!
# If you get this error:
# "ERROR multiple chr/scaf detected, VCF should contain a single chr/scaf for plotting Execution halted"
# but you are sure that every file is chromosome specific, check the vcf as it could just be a chromosome with <=1 
# site designation and therefore there isn't enough info in the vcf (since we subsetted with site des)! No worries if this is the case :)

# Install & load necessary packages ####

# List of packages needed
packages <- c("optparse", # v 1.7.5, allow use of command line flags
              "remotes", # v 2.5.0
              "vcfR",# v 1.15.0, read.vcfR
              "GenotypePlot" # v 0.2.1, 
) 


installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  
  to_install <- packages[!installed_packages] # necessary packages that aren't installed yet
  
  special_packages <- c("GenotypePlot") # this package need to be installed specially!
  
  install.packages(setdiff(to_install, special_packages), repos = "https://cloud.r-project.org/") # setdiff pulls package names that aren't in special_packages vector but are in to_install
  # won't run in command line without setting repos, CRAN mirror error message!
  
  # install this pacakge differently!
  if ("GenotypePlot" %in% to_install) { # if genotype plot needs to be installed
    remotes::install_github("JimWhiting91/genotype_plot") # then do it this way
  }
}

# Load the packages
invisible(lapply(packages, library, character.only = TRUE))


# Run from command line! ####

flag_list <- list(
  make_option(c("--VCF"), type = "character", default = NULL, 
              help = "Path to vcf!", metavar = "file_path"),
  
  make_option(c("--width"), type = "integer", default = 10, 
              help = "Width of output pdf [default= %default]", metavar = "integer"),
  
  make_option(c("--height"), type = "integer", default = 8, 
              help = "Height of output pdf [default= %default]", metavar = "integer")
)

# Set the flag list
opt_parser <- OptionParser(option_list = flag_list)
opt <- parse_args(opt_parser)

# Imports ####

# call the input file "vcf"
vcf <- opt$VCF

# Set the popmap ####
# popmap = two column data frame with column 1 for individual IDs as they appear in the VCF and column 2 for pop labels
our_popmap <- data.frame(ind = c("K006", "K010", "K011", "K015", "K019", "K021", "K022", "K026", "K029", "K031", "K032", "034", "K035", "K039", "K071", "K073", "K083", "K085", "K089", "K090", "K096", 
                                 "C060", "C065", "C067", "C068", "C080", "C081", "C088", "C090", "C094",
                                 "8L19786", "8L64869", "8L89915", "8M31651", "8M71932", "8M72455", "8N05890", "8N06612", "8N73248", "8N73604", 
                                 "M010", "M012", "M019", "M022", "M028", "M030", "M031", "M033", "M036", "M037",
                                 "8934547", "8L19766", "8L52141", "8L52830", "8N05240",
                                 "N33", "PM001-Malta", "PM002-Malta", "PM003", "PM007", "PM2", "PM4", "PM5", "PM6", 
                                 "S029", "S040", "S045", "S046", "S049", "S052", "S058", "S059", "S064", "S065", 
                                 "Lesina_280", "Lesina_281", "Lesina_282", "Lesina_285", "Lesina_286", "Lesina_287", "Lesina_288", "Lesina_289", "Lesina_292", "Lesina_295"),
                         pop = c(rep("Corsica", 21), rep("Crete", 9), rep("House", 10), rep("Malta", 10), rep("House", 5), rep("Tree", 9), rep("Sicily", 10), rep("Spanish", 10)),
                         stringsAsFactors = FALSE)



# Generate plot! ####
# Using a loop to loop through VCFs in directory doesn't work well in R! Outputs empty pdfs.
# Instead, use R script for one vcf at a time but creat a bash loop to run for each vcf

cat("\nProcessing vcf:", vcf, "\n\n") # Print to user screen to see which vcf is being processed (helpful if in bash loop)
  
my_vcf <- read.vcfR(vcf) # Read the vcf file in 

# pulled straight from github readme
new_plot <- genotype_plot(vcf_object  = my_vcf,
                          popmap = our_popmap,                              
                          cluster        = FALSE,                           
                          snp_label_size = 10000,                          
                          colour_scheme=c("#d4b9da","#e7298a","#980043")) 


vcf_basename <- tools::file_path_sans_ext(basename(vcf)) # pull the basename, e.g CM071427.1_subset_cutoff_0.99
output_file <- paste0(vcf_basename, "_genotypePlot.pdf") # CM071427.1_subset_cutoff_0.99_genotypePlot.pdf

pdf(output_file, width = opt$width, height = opt$height) # name the pdf!

combine_genotype_plot(new_plot) # plot!

dev.off()




