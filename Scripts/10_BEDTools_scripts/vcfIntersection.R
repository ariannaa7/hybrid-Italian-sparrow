# Load packages ####
library(tidyverse) # data frame manipulation

# Import Data ####

# in 0.90 250,695 sites total and in 0.99 137,411 sites total
# VCFs subsetted based on designation files -> VCFs intersected with gff file to produce subsetted vcfs based on feature type
# The counted the number of vcf sites for each cutoff and feature e.g. cat exon_sites/*vcf | grep -v "#" | wc -l
# There will be overlap between genes & exons! And between intergenic sites and sites within 10,000

  VCF_intersection_stats = data.frame(
    cutoff = c(0.90,0.99),
    exon_sites_count = c(15221,9176),
    gene_sites_count = c(142079,81780),
    intergenic_sites_count = c(108616,55631),
    nearGene_sites_count = c(41419,23992)
  )
  
  VCF_intersection_stats_autosomes = data.frame(
    cutoff = c(0.90,0.99),
    exon_sites_count = c(5676,2271),
    gene_sites_count = c(50898,17473),
    intergenic_sites_count = c(31641,7600),
    nearGene_sites_count = c(11261,3525)
  )
  
  VCF_intersection_stats_Z = data.frame(
    cutoff = c(0.90,0.99),
    exon_sites_count = c(9545,6905),
    gene_sites_count = c(91181,64307),
    intergenic_sites_count = c(76975,48031),
    nearGene_sites_count = c(30158,20467)
  )
  
  # Long format
  VCF_intersection_stats_long_format <- function(stats_df) {
    long_df = stats_df %>%
      pivot_longer(cols = -cutoff, names_to = "count_type", values_to = "count")
    
    return(long_df)
  }
  
  VCF_intersection_stats_long <- VCF_intersection_stats_long_format(VCF_intersection_stats)
  VCF_intersection_stats_autosomes_long <- VCF_intersection_stats_long_format(VCF_intersection_stats_autosomes)
  VCF_intersection_stats_Z_long <- VCF_intersection_stats_long_format(VCF_intersection_stats_Z)
  

# Bar plot function
  plot_intersectionStats <- function(VCF_intersection_stats_long_df, title) { # title should be string with ""
    plot <-  ggplot(VCF_intersection_stats_long_df, aes(x = as.factor(cutoff), y = count, fill = count_type)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = "Allele frequency cutoffs", y = "Number of designated sites falling within this region", fill = "Feature Type") + 
      scale_y_continuous(labels = scales::comma_format()) +
      scale_fill_manual(values = c("exon_sites_count" = "seagreen", "gene_sites_count" = "sienna", "intergenic_sites_count" = "slateblue", "nearGene_sites_count" = "palevioletred"),
                        labels = c("exon_sites_count" = "Exons", "gene_sites_count" = "Genes", "intergenic_sites_count" = "Intergenic Sites", "nearGene_sites_count" = "Sites within 10,000 bps of a gene")) +
      ggtitle(title) +
      scale_x_discrete(labels = c("0.9" = "0.90")) +
      theme_minimal() +
      theme(legend.title = element_text(hjust = 0.5)) # center the legend title
    
    return(plot)
  } 
  
  # Plot!
  plot_intersectionStats(VCF_intersection_stats_long, "Site designation counts per cutoff")
  plot_intersectionStats(VCF_intersection_stats_autosomes_long, "Site designation counts per cutoff - autosomes")
  plot_intersectionStats(VCF_intersection_stats_Z_long, "Site designation counts per cutoff - Z")
