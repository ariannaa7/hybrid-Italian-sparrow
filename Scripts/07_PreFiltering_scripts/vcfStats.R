# The following comes from the https://speciationgenomics.github.io/filtering_vcfs/ guide to filtering VCFs, with minor adjustments

# Load packages -----
library(tidyverse)

# Set Working directory ----
#setwd("~/BINP37/VCF_filteringThreshold_stats/VCFtools")

# sparrow_biallelics_subset ----
  # Variant quality ----
  var_qual <- read_delim("sparrow_biallelics_subset.lqual", delim = "\t",
                         col_names = c("chr", "pos", "qual"), skip = 1)
  
  var_qual_plot <- ggplot(var_qual, aes(qual)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light() 
  
  plot(var_qual_plot)

  # Variant mean depth ----
  var_depth <- read_delim("sparrow_biallelics_subset.ldepth.mean", delim = "\t",
                          col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
  
  var_depth_plot <- ggplot(var_depth, aes(mean_depth)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(var_depth_plot)
  
  # Pull summary stats for the mean depth column
  summary(var_depth$mean_depth) # most variants have a depth of 8.49-10.01x (1st, med, mean, 3rd)
  
  # Change x range for better visualization
  var_depth_plot + theme_light() + xlim(0, 100)
  
  # Variant missingness ----
  var_miss <- read_delim("sparrow_biallelics_subset.lmiss", delim = "\t",
                         col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  var_miss_plot <- ggplot(var_miss, aes(fmiss)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_miss_plot)
  
  # Pull summary stats for fraction missingness per site?
  summary(var_miss$fmiss) # most sites have minimal missing data, close to none
  
  # Minor allele frequency ----
  var_freq <- read_delim("sparrow_biallelics_subset.frq", delim = "\t",
                         col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  
  # calculate MAF at each site
  var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
  
  var_freq_plot <- ggplot(var_freq, aes(maf)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_freq_plot)
  
  # Pull summary stats for MAF per site
  summary(var_freq$maf)
  
  # Mean depth per individual ----
  ind_depth <- read_delim("sparrow_biallelics_subset.idepth", delim = "\t",
                          col_names = c("ind", "nsites", "depth"), skip = 1)
  
  ind_depth_plot <- ggplot(ind_depth, aes(depth)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_depth_plot)
  
  # Proportion missing data per individual ----
  ind_miss <- read_delim("sparrow_biallelics_subset.imiss", delim = "\t",
                          col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  ind_miss_plot <- ggplot(ind_miss, aes(fmiss)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_miss_plot)

  # Heterozygosity and inbreeding coefficient per individual ----
  ind_het <- read_delim("sparrow_biallelics_subset.het", delim = "\t",
                        col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)

  ind_het_plot <- ggplot(ind_het, aes(f)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_het_plot)


  
# sparrow_biallelics_subset_AUTOSOMES ----
  
  # Variant quality ----
  var_qual <- read_delim("sparrow_biallelics_subset_AUTOSOMES.lqual", delim = "\t",
                        col_names = c("chr", "pos", "qual"), skip = 1)
  
  var_qual_plot <- ggplot(var_qual, aes(qual)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light() 
  
  plot(var_qual_plot)
  
  var_qual_plot + ggtitle("Variant Quality - Autosomes")
  var_qual_plot + xlim(0,10000) + ggtitle("Variant Quality - Autosomes, smaller xlim")
  
  
  # Variant mean depth ----
  var_depth <- read_delim("sparrow_biallelics_subset_AUTOSOMES.ldepth.mean", delim = "\t",
                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
  
  var_depth_plot <- ggplot(var_depth, aes(mean_depth)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(var_depth_plot)
  
  # Pull summary stats for the mean depth column
  summary(var_depth$mean_depth) # most variants have a depth of 8.52-10.03x (1st, med, mean, 3rd)
  
  # Change x range for better visualization
  var_depth_plot + theme_light() + xlim(0, 15) + 
    ggtitle(label = "Variant Mean Depth - Autosomes",
            subtitle = "min = 0.02, 1stQ = 8.52, Med = 9.45, Mean = 9.09, 3rdQ = 10.03, Max = 152.5")

  # Variant missingness ----
  var_miss <- read_delim("sparrow_biallelics_subset_AUTOSOMES.lmiss", delim = "\t",
                        col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  var_miss_plot <- ggplot(var_miss, aes(fmiss)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_miss_plot) + xlim(0,0.2) +
    ggtitle(label = "Variant Missingness - Autosomes",
            subtitle = "min = 0.00, 1stQ = 0.00, Med = 0.00, Mean = 0.02, 3rdQ = 0.00, Max = 0.98")
  
  
  # Pull summary stats for fraction missingness per site?
  summary(var_miss$fmiss) # most sites have minimal missing data, close to none
  
  # Minor allele frequency ----
  var_freq <- read_delim("sparrow_biallelics_subset_AUTOSOMES.frq", delim = "\t",
                        col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  
  # calculate MAF at each site
  var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
  
  var_freq_plot <- ggplot(var_freq, aes(maf)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_freq_plot) +
    ggtitle(label = "MAF - Autosomes",
            subtitle = "min = 0.00, 1stQ = 0.008, Med = 0.03, Mean = 0.08, 3rdQ = 0.08")
  
  
  # Pull summary stats for MAF per site
  summary(var_freq$maf)
  
  # Mean depth per individual ----
  ind_depth <- read_delim("sparrow_biallelics_subset_AUTOSOMES.idepth", delim = "\t",
                         col_names = c("ind", "nsites", "depth"), skip = 1)
  
  ind_depth_plot <- ggplot(ind_depth, aes(depth)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_depth_plot) +
    ggtitle(label = "Individual Mean Depth - Autosomes")
  
  summary(ind_depth$depth)
  
  # Proportion missing data per individual ----
  ind_miss <- read_delim("sparrow_biallelics_subset_AUTOSOMES.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  ind_miss_plot <- ggplot(ind_miss, aes(fmiss)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_miss_plot) +
    ggtitle(label = "Individual Missingness - Autosomes")
  
  # Heterozygosity and inbreeding coefficient per individual ----
  ind_het <- read_delim("sparrow_biallelics_subset_AUTOSOMES.het", delim = "\t",
                       col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
  
  ind_het_plot <- ggplot(ind_het, aes(f)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_het_plot)
  
  
  
# sparrow_biallelics_subset_SEX ----
  
  # Variant quality ----
  var_qual <- read_delim("sparrow_biallelics_subset_SEX.lqual", delim = "\t",
                        col_names = c("chr", "pos", "qual"), skip = 1)
  
  var_qual_plot <- ggplot(var_qual, aes(qual)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light() 
  
  plot(var_qual_plot)
  
  # Variant mean depth ----
  var_depth <- read_delim("sparrow_biallelics_subset_SEX.ldepth.mean", delim = "\t",
                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
  
  var_depth_plot <- ggplot(var_depth, aes(mean_depth)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(var_depth_plot)
  
  # Pull summary stats for the mean depth column
  summary(var_depth$mean_depth) # most variants have a depth of 7.84-9.75x (1st, med, mean, 3rd)
  
  # Change x range for better visualization
  var_depth_plot + theme_light() + xlim(0, 20) +
    ggtitle(label = "Variant Mean Depth - Z & W chromosomes",
            subtitle = "1stQ = 7.85, Med = 9.18, Mean = 8.12, 3rdQ = 9.75")
  
  
  # Variant missingness ----
  var_miss <- read_delim("sparrow_biallelics_subset_SEX.lmiss", delim = "\t",
                        col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  var_miss_plot <- ggplot(var_miss, aes(fmiss)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_miss_plot)
  
  # Pull summary stats for fraction missingness per site?
  summary(var_miss$fmiss) # slightly more missingness that autosome, still minimal
  
  # Minor allele frequency ----
  var_freq <- read_delim("sparrow_biallelics_subset_SEX.frq", delim = "\t",
                        col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  
  # calculate MAF at each site
  var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
  
  var_freq_plot <- ggplot(var_freq, aes(maf)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_freq_plot)
  
  # Pull summary stats for MAF per site
  summary(var_freq$maf)
  
  # Mean depth per individual ----
  ind_depth <- read_delim("sparrow_biallelics_subset_SEX.idepth", delim = "\t",
                         col_names = c("ind", "nsites", "depth"), skip = 1)
  
  ind_depth_plot <- ggplot(ind_depth, aes(depth)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_depth_plot)
  
  # Proportion missing data per individual ----
  ind_miss <- read_delim("sparrow_biallelics_subset_SEX.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  ind_miss_plot <- ggplot(ind_miss, aes(fmiss)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_miss_plot)
  
  # Heterozygosity and inbreeding coefficient per individual ----
  ind_het <- read_delim("sparrow_biallelics_subset_SEX.het", delim = "\t",
                       col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
  
  ind_het_plot <- ggplot(ind_het, aes(f)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_het_plot)
  
  
  
# sparrow_biallelics_subset_wChrom ----
  # Variant quality ----
  var_qual <- read_delim("sparrow_biallelics_subset_wChrom.lqual", delim = "\t",
                        col_names = c("chr", "pos", "qual"), skip = 1)
  
  var_qual_plot <- ggplot(var_qual, aes(qual)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light() 
  
  plot(var_qual_plot)
  
  # Variant mean depth ----
  var_depth <- read_delim("sparrow_biallelics_subset_wChrom.ldepth.mean", delim = "\t",
                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
  
  var_depth_plot <- ggplot(var_depth, aes(mean_depth)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(var_depth_plot)
  
  # Pull summary stats for the mean depth column
  summary(var_depth$mean_depth) # most variants have a depth of 0.2-1.05x (1st, med, mean, 3rd)
  
  # Change x range for better visualization
  var_depth_plot + theme_light() + xlim(0, 100)
  
  # Variant missingness ----
  var_miss <- read_delim("sparrow_biallelics_subset_wChrom.lmiss", delim = "\t",
                        col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  var_miss_plot <- ggplot(var_miss, aes(fmiss)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_miss_plot)
  
  # Pull summary stats for fraction missingness per site?
  summary(var_miss$fmiss) # more missingness that autosome
  
  # Minor allele frequency ----
  var_freq <- read_delim("sparrow_biallelics_subset_wChrom.frq", delim = "\t",
                        col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  
  # calculate MAF at each site
  var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
  
  var_freq_plot <- ggplot(var_freq, aes(maf)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_freq_plot)
  
  # Pull summary stats for MAF per site
  summary(var_freq$maf)
  
  # Mean depth per individual ----
  ind_depth <- read_delim("sparrow_biallelics_subset_wChrom.idepth", delim = "\t",
                         col_names = c("ind", "nsites", "depth"), skip = 1)
  
  ind_depth_plot <- ggplot(ind_depth, aes(depth)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_depth_plot)
  
  # Proportion missing data per individual ----
  ind_miss <- read_delim("sparrow_biallelics_subset_wChrom.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  ind_miss_plot <- ggplot(ind_miss, aes(fmiss)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_miss_plot)
  
  # Heterozygosity and inbreeding coefficient per individual ----
  ind_het <- read_delim("sparrow_biallelics_subset_wChrom.het", delim = "\t",
                       col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
  
  ind_het_plot <- ggplot(ind_het, aes(f)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_het_plot)
  


# sparrow_biallelics_subset_zChrom ----
  
  # Variant quality ----
  var_qual <- read_delim("sparrow_biallelics_subset_zChrom.lqual", delim = "\t",
                        col_names = c("chr", "pos", "qual"), skip = 1)
  
  var_qual_plot <- ggplot(var_qual, aes(qual)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light() 
  
  plot(var_qual_plot)
  
  # Variant mean depth ----
  var_depth <- read_delim("sparrow_biallelics_subset_zChrom.ldepth.mean", delim = "\t",
                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
  
  var_depth_plot <- ggplot(var_depth, aes(mean_depth)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(var_depth_plot)
  
  # Pull summary stats for the mean depth column
  summary(var_depth$mean_depth) # most variants have a depth of 8.69-9.82x (1st, med, mean, 3rd)
  
  # Change x range for better visualization
  var_depth_plot + theme_light() + xlim(0, 100)
  
  # Variant missingness ----
  var_miss <- read_delim("sparrow_biallelics_subset_zChrom.lmiss", delim = "\t",
                        col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  var_miss_plot <- ggplot(var_miss, aes(fmiss)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_miss_plot)
  
  # Pull summary stats for fraction missingness per site?
  summary(var_miss$fmiss) # minimal missingness
  
  # Minor allele frequency ----
  var_freq <- read_delim("sparrow_biallelics_subset_zChrom.frq", delim = "\t",
                        col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  
  # calculate MAF at each site
  var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
  
  var_freq_plot <- ggplot(var_freq, aes(maf)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_freq_plot)
  
  # Pull summary stats for MAF per site
  summary(var_freq$maf)
  
  # Mean depth per individual ----
  ind_depth <- read_delim("sparrow_biallelics_subset_zChrom.idepth", delim = "\t",
                         col_names = c("ind", "nsites", "depth"), skip = 1)
  
  ind_depth_plot <- ggplot(ind_depth, aes(depth)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_depth_plot)
  
  # Proportion missing data per individual ----
  ind_miss <- read_delim("sparrow_biallelics_subset_zChrom.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  ind_miss_plot <- ggplot(ind_miss, aes(fmiss)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_miss_plot)
  
  # Heterozygosity and inbreeding coefficient per individual ----
  ind_het <- read_delim("sparrow_biallelics_subset_zChrom.het", delim = "\t",
                       col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
  
  ind_het_plot <- ggplot(ind_het, aes(f)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_het_plot)
  
  
  
# sparrow_biallelics_subset_wChrom_REDO ----
  # Variant quality ----
  var_qual <- read_delim("sparrow_biallelics_subset_wChrom_REDO.lqual", delim = "\t",
                        col_names = c("chr", "pos", "qual"), skip = 1)
  
  var_qual_plot <- ggplot(var_qual, aes(qual)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light() 
  
  plot(var_qual_plot) + 
    xlim(0,10000) + ggtitle("Variant Quality - W chromosome, smaller xlim")
  
  # Variant mean depth ----
  var_depth <- read_delim("sparrow_biallelics_subset_wChrom_REDO.ldepth.mean", delim = "\t",
                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
  
  var_depth_plot <- ggplot(var_depth, aes(mean_depth)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(var_depth_plot) 
  
  # Pull summary stats for the mean depth column
  summary(var_depth$mean_depth) # most variants have a depth of 0.2-1.07x (1st, med, mean, 3rd)
  
  # Change x range for better visualization
  var_depth_plot + theme_light() + xlim(0, 5) +
    ggtitle(label = "Variant Mean Depth - W chromosome",
            subtitle = "min = 0.008, 1stQ = 0.22, Med = 0.354, Mean = 1.25, 3rdQ = 1.07, Max = 145.11")
  
  
  
  # Variant missingness ----
  var_miss <- read_delim("sparrow_biallelics_subset_wChrom_REDO.lmiss", delim = "\t",
                        col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  var_miss_plot <- ggplot(var_miss, aes(fmiss)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_miss_plot) +
    ggtitle(label = "Variant Missingness - W chromosome",
            subtitle = "min = 0.00, 1stQ = 0.56, Med = 0.90, Mean = 0.73, 3rdQ = 0.96, Max = 0.99")
  
  
  # Pull summary stats for fraction missingness per site?
  summary(var_miss$fmiss) # more missingness that autosome
  
  # Minor allele frequency ----
  var_freq <- read_delim("sparrow_biallelics_subset_wChrom_REDO.frq", delim = "\t",
                        col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  
  # calculate MAF at each site
  var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
  
  var_freq_plot <- ggplot(var_freq, aes(maf)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_freq_plot) +
    ggtitle(label = "MAF - W chromosome",
            subtitle = "min = 0.00, 1stQ = 0.00, Med = 0.08, Mean = 0.12, 3rdQ = 0.2")
  
  
  # Pull summary stats for MAF per site
  summary(var_freq$maf)
  
  # Mean depth per individual ----
  ind_depth <- read_delim("sparrow_biallelics_subset_wChrom_REDO.idepth", delim = "\t",
                         col_names = c("ind", "nsites", "depth"), skip = 1)
  
  ind_depth_plot <- ggplot(ind_depth, aes(depth)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_depth_plot) +
    ggtitle(label = "Individual Mean Depth - W")
  
  # Proportion missing data per individual ----
  ind_miss <- read_delim("sparrow_biallelics_subset_wChrom_REDO.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  ind_miss_plot <- ggplot(ind_miss, aes(fmiss)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_miss_plot) +
    ggtitle(label = "Individual Missingness - W")
  
  # Heterozygosity and inbreeding coefficient per individual ----
  ind_het <- read_delim("sparrow_biallelics_subset_wChrom_REDO.het", delim = "\t",
                       col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
  
  ind_het_plot <- ggplot(ind_het, aes(f)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_het_plot)
  
  
  
# sparrow_biallelics_subset_zChrom_REDO ----
  # Variant quality ----
  var_qual <- read_delim("sparrow_biallelics_subset_zChrom_REDO.lqual", delim = "\t",
                        col_names = c("chr", "pos", "qual"), skip = 1)
  
  var_qual_plot <- ggplot(var_qual, aes(qual)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light() 
  
  plot(var_qual_plot) + xlim(0,10000) +
    ggtitle("Variant Quality - Z chromosome, smaller xlim")
  
  # Variant mean depth ----
  var_depth <- read_delim("sparrow_biallelics_subset_zChrom_REDO.ldepth.mean", delim = "\t",
                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
  
  var_depth_plot <- ggplot(var_depth, aes(mean_depth)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(var_depth_plot)
  
  # Pull summary stats for the mean depth column
  summary(var_depth$mean_depth) # most variants have a depth of 8.6-9.8x (1st, med, mean, 3rd)
  
  # Change x range for better visualization
  var_depth_plot + theme_light() + xlim(0, 13) +
    ggtitle(label <- "Variant Mean Depth - Z chromosome",
            subtitle = "min = 0.02, 1stQ = 8.68, Med = 9.34, Mean = 9.30, 3rdQ = 9.84, Max = 207.19")
  
  var_depth_plot + theme_light() + xlim(0, 12)
  
  # Variant missingness ----
  var_miss <- read_delim("sparrow_biallelics_subset_zChrom_REDO.lmiss", delim = "\t",
                        col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  var_miss_plot <- ggplot(var_miss, aes(fmiss)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_miss_plot) +
    ggtitle(label = "Variant Missingness - Z chromosome",
            subtitle = "min = 0.00, 1stQ = 0.00, Med = 0.00, Mean = 0.01, 3rdQ = 0.00, Max = 0.97")
  
  
  # Pull summary stats for fraction missingness per site?
  summary(var_miss$fmiss) # more missingness that autosome
  
  # Minor allele frequency ----
  var_freq <- read_delim("sparrow_biallelics_subset_zChrom_REDO.frq", delim = "\t",
                        col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  
  # calculate MAF at each site
  var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
  
  var_freq_plot <- ggplot(var_freq, aes(maf)) + 
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
    theme_light()
  
  plot(var_freq_plot) +
    ggtitle(label = "MAF - Z chromosome",
            subtitle = "min = 0.00, 1stQ = 0.008, Med = 0.04, Mean = 0.09, 3rdQ = 0.09")
  
  
  # Pull summary stats for MAF per site
  summary(var_freq$maf)
  
  # Mean depth per individual ----
  ind_depth <- read_delim("sparrow_biallelics_subset_zChrom_REDO.idepth", delim = "\t",
                         col_names = c("ind", "nsites", "depth"), skip = 1)
  
  ind_depth_plot <- ggplot(ind_depth, aes(depth)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_depth_plot) +
    ggtitle(label = "Individual Mean Depth - Z")
  
  # Proportion missing data per individual ----
  ind_miss <- read_delim("sparrow_biallelics_subset_zChrom_REDO.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  ind_miss_plot <- ggplot(ind_miss, aes(fmiss)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_miss_plot) + ggtitle(label = "Individual Missingness - Z")
  
   
  # Heterozygosity and inbreeding coefficient per individual ----
  ind_het <- read_delim("sparrow_biallelics_subset_zChrom_REDO.het", delim = "\t",
                       col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
  
  ind_het_plot <- ggplot(ind_het, aes(f)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light()
  
  plot(ind_het_plot)
  
