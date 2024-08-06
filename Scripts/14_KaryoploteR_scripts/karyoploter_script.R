# The following code is used to visual how derived sites are distributed across chromosomes in the genome
# Has been adapted from code provided at https://bernatgel.github.io/karyoploter_tutorial/

# Load Packages ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("karyoploteR")
library(karyoploteR)
library(dplyr)
library(tidyr)
library(magrittr)
library(rtracklayer)

# GFF file, download gff file (refseq) from https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_036417665.1/ ####

gff.file <- "genomic.gff" # create a variable that stores the path to the gff file
lines <- readLines(gff.file) # read in the gff line by line

#The lines with the standard chromosomes start with "##sequence-region NC_0", read them in
# Scaffolds and super chroms start start with "NW", unreadable plot, leave out, not many sites there anyway
chrom_lines <- lines[grepl(lines, pattern = "##sequence-region NC_0")] # store "##sequence-region NC_087474.1 1 161085800""

# Split by space " " and turn these gff lines into a dataframe
gff_df <- data.frame(do.call(rbind, strsplit(chrom_lines, split = " ")))
gff_df[,3] <- as.numeric(as.character(gff_df[,3])) # stand information as numeric
gff_df[,4] <- as.numeric(as.character(gff_df[,4])) # chrom length/size as numeric

#pdom.genome <- toGRanges(gff_df[,c(1,2,3)]) # converts gff to GRanges (genomic ranges) object, P. domesticus genome

# Rename the columns of the gff to prepare for merging with seq_report
gff_df <- gff_df %>% rename(headerLine = X1, # unecessary column!
                            RefSeq = X2, 
                            strand = X3, 
                            Size_gff = X4) 

# Sequence report file, download from https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_036417665.1/ ####

# We need the seq report because the gff is using RefSeq accession, but all of our data is based on GenBank accession
# Use the sequence report to change refseq to genbank
seq_report <- read.csv("sequence_report.tsv", header = T, sep = "\t") %>% 
  select("Chromosome.name", "GenBank.seq.accession","RefSeq.seq.accession", "Seq.length") %>%
  rename( Chromosome = Chromosome.name,
          GenBank = GenBank.seq.accession,
          RefSeq = RefSeq.seq.accession,
          Size = Seq.length # name size differently than size from gff, this size is integer and will be treated differently (lengths are exact same between gff & seq report)
  )


# Now, join the gff_df with the seq_report df to change refseq accessions to genbank accessions
gff_df <- full_join(gff_df, seq_report) # joins by refSeq accession
gff_df <- gff_df %>% 
  select(GenBank, strand, Size, Chromosome) %>%
  filter(GenBank != "KM078784.1") %>% # Exclude mitochondrial genome!
  filter(GenBank != "CM071464.1") %>% # Exclude W chromosome !
  drop_na()

pdom.genome <- toGRanges(gff_df[,c(1,2,3)]) # convert to a genomic range object (chrom, strand, and size)

# Import BED files! ####
SD_noPrune <- import("trimmed_SD_nonPruned_allChroms.bed") # Trim bed files to only include first 4 cols, cat *.bed | cut -f 1-4 >
HD_noPrune <- import("trimmed_HD_nonPruned_allChroms.bed") # Issues reading file in if you keep all cols, the other ones are unecessary for this analysis

SD_50 <- import("trimmed_SD_window50_allChroms.bed") # Trim bed files to only include first 4 cols, cat *.bed | cut -f 1-4 >
HD_50 <- import("trimmed_HD_window50_allChroms.bed") # Issues reading file in if you keep all cols, the other ones are unecessary for this analysis

SD_200 <- import("trimmed_SD_window200_allChroms.bed") # Trim bed files to only include first 4 cols, cat *.bed | cut -f 1-4 >
HD_200 <- import("trimmed_HD_window200_allChroms.bed") # Issues reading file in if you keep all cols, the other ones are unecessary for this analysis

SD_450 <- import("trimmed_SD_window450_allChroms.bed") # Trim bed files to only include first 4 cols, cat *.bed | cut -f 1-4 >
HD_450 <- import("trimmed_HD_window450_allChroms.bed") # Issues reading file in if you keep all cols, the other ones are unecessary for this analysis

# The bed files are already loaded in as genomic range objects!!

# Plot

# No pruning
png(filename = "derSites_acrossChrom_noPruning.png", width = 3500, height = 2000, res = 300)

kp_noPrun <- plotKaryotype(genome=pdom.genome, plot.type = 2) # Horizontal ideograms with two data panels, one above and one below them
kpPlotRegions(kp_noPrun, data=SD_noPrune, avoid.overlapping = FALSE, data.panel=1, col="#BE441E") # SD above line
kpPlotRegions(kp_noPrun, data=HD_noPrune, avoid.overlapping = FALSE, data.panel=2, col = "#4B68D9") # HD below line
legend("right", fill = c("#BE441E", "#4B68D9"), legend = c("Spanish Derived", "House Derived"), cex = 1.1, bty="n")

dev.off()

# 50 kb window
png(filename = "derSites_acrossChrom_50kbWindow.png", width = 3500, height = 2000, res = 300)

kp_50 <- plotKaryotype(genome=pdom.genome, plot.type = 2) # Horizontal ideograms with two data panels, one above and one below them
kpPlotRegions(kp_50, data=SD_50, avoid.overlapping = FALSE, data.panel=1, col="#BE441E") # SD above line
kpPlotRegions(kp_50, data=HD_50, avoid.overlapping = FALSE, data.panel=2, col = "#4B68D9") # HD below line
legend("right", fill = c("#BE441E", "#4B68D9"), legend = c("Spanish Derived", "House Derived"), cex = 1.1, bty="n")

dev.off()

# 200 kb window
png(filename = "derSites_acrossChrom_200kbWindow.png", width = 3500, height = 2000, res = 300)

kp_200 <- plotKaryotype(genome=pdom.genome, plot.type = 2) # Horizontal ideograms with two data panels, one above and one below them
kpPlotRegions(kp_200, data=SD_200, avoid.overlapping = FALSE, data.panel=1, col="#BE441E") # SD above line
kpPlotRegions(kp_200, data=HD_200, avoid.overlapping = FALSE, data.panel=2, col = "#4B68D9") # HD below line
legend("right", fill = c("#BE441E", "#4B68D9"), legend = c("Spanish Derived", "House Derived"), cex = 1.1, bty="n")

dev.off()

# 450 kb window
png(filename = "derSites_acrossChrom_450kbWindow.png", width = 3500, height = 2000, res = 300)

kp_450 <- plotKaryotype(genome=pdom.genome, plot.type = 2) # Horizontal ideograms with two data panels, one above and one below them
kpPlotRegions(kp_450, data=SD_450, avoid.overlapping = FALSE, data.panel=1, col="#BE441E") # SD above line
kpPlotRegions(kp_450, data=HD_450, avoid.overlapping = FALSE, data.panel=2, col = "#4B68D9") # HD below line
legend("right", fill = c("#BE441E", "#4B68D9"), legend = c("Spanish Derived", "House Derived"), cex = 1.1, bty="n")

dev.off()
