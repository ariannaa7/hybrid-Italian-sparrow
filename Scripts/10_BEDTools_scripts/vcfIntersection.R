# Load packages ####
library(tidyverse) # data frame manipulation

# Import Data ####

# in 0.90 250,695 sites total and in 0.99 137,411 sites total
# VCFs subsetted based on designation files -> VCFs intersected with gff file to produce subsetted vcfs based on feature type
# The counted the number of vcf sites for each cutoff and feature e.g. cat exon_sites/*vcf | grep -v "#" | wc -l

  GFF_intersection_stats = data.frame(
    cutoff = c(0.90,0.99),
    exon_sites_count = c(15221,9176),
    gene_sites_count = c(142079,81780),
    intergenic_sites_count = c(108616,55631),
    nearGene_sites_count = c(41419,23992)
  )
  
  
  # Long format
  GFF_intersection_stats_long = GFF_intersection_stats %>%
    pivot_longer(cols = -cutoff, names_to = "count_type", values_to = "count")

# Stacked bar plot
  ggplot(GFF_intersection_stats_long, aes(x = as.factor(cutoff), y = count, fill = count_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Allele frequency cutoffs", y = "Number of designated sites falling within this region", fill = "Feature Type") + 
    scale_y_continuous(labels = scales::comma_format()) +
    scale_fill_manual(values = c("exon_sites_count" = "seagreen", "gene_sites_count" = "sienna", "intergenic_sites_count" = "slateblue", "nearGene_sites_count" = "palevioletred"),
                      labels = c("exon_sites_count" = "Exons", "gene_sites_count" = "Genes", "intergenic_sites_count" = "Intergenic Sites", "nearGene_sites_count" = "Sites within 10,000 bps of a gene")) +
    #ggtitle("Site designation counts per cutoff") +
    scale_x_discrete(labels = c("0.9" = "0.90")) +
    theme_minimal() +
    theme(legend.title = element_text(hjust = 0.5)) # center the legend title
  