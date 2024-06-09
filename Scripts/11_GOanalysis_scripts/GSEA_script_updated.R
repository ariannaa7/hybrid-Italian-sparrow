#####
# Script was created/used by Kalle & Rachel
# Updates made by Arianna: 
# to intall topGO and Rgraphviz, one must have BiocManager installed!
# rm(list = objects())
# Removed unecessary comment to add header

# This is a  script to take in two files that have been pre processed
# 1) you need an annotation file that has a list of all genes tab then 
# comma separated GO terms
# 2) you need a list of the candidate genes from that list you want to assess 
# for being enriched compared to the genome.
#
# 1st argument is the annotation file (header: GeneID\tGO_ID)
# 2nd argument is the candidate set list (header: geneID)
#
# note that the current setting for GO node size is 5, which 
# is hard coded below for more robust results.
#
# example run
# Rscript GSEA_run_script.R annotationFile.tsv candidateSet.tsv

# descriptions
# http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html

# new general run

rm(list = objects()) #clears all variables
graphics.off() #close all figures
packages <- c("BiocManager","topGO","Rgraphviz", "openxlsx") # BiocManager is needed to install topGO and Rgraphviz

# BiocManager v 1.30.23
# openxlsx v 4.2.5.2


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  
  to_install <- packages[!installed_packages] # necessary packages that aren't installed yet
  
  bioc_packages <- c("topGO", "Rgraphviz") # these packages need to be installed specially!
  
  install.packages(setdiff(to_install, bioc_packages), repos = "https://cloud.r-project.org/") # setdiff pulls package names that aren't in bioc_packages vector but are in to_install
  # so, will install BiocManager and openxlsx if not installed already
  # won't run in command line without setting repos, CRAN mirror error message!
  
  # install these pacakges differently!
  BiocManager::install(intersect(to_install, bioc_packages)) # set diff pulls package names that are in bioc_packages AND to_install
  
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
##########
# set variables
# making general script 
# Rscript myscript.R batch.csv
# and invoke these in the myscript.R
args <- commandArgs(TRUE)
# dataset <- read.table(args[1],header=FALSE,sep=",",skip=1)

annotations=args[1]
candidate_list=args[2]
# here I will be only analyzing GO terms with at least 5 members,
# as this yield more stable results.
node_size=5

#### CORE BP #####
GO_category="BP" # biological process
geneID2GO <- readMappings(file = annotations)  
geneUniverse <- names(geneID2GO) 

genesOfInterest.bv <- read.table(candidate_list,header=TRUE)

genesOfInterest.bv <- as.character(genesOfInterest.bv$geneID) 
geneList.bv <- factor(as.integer(geneUniverse %in% genesOfInterest.bv))
names(geneList.bv) <- geneUniverse

myGOdata.bv <- new("topGOdata", description="Candidate genes", ontology=GO_category, allGenes=geneList.bv,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = node_size)

# STATS#
# each GO term is tested independently, not taking the GO hierarchy into account
resultClassic <- runTest(myGOdata.bv, algorithm="classic", statistic="fisher")
# elim method processes the GO terms by traversing the GO hierarchy from bottom to top, 
# ie. it first assesses the most specific (bottom-most) GO terms, and proceeds later 
# to more general (higher) GO terms. When it assesses a higher (more general) GO term, 
# it discards any genes that are annotated with significantly enriched descendant 
# GO terms (considered significant using a pre-defined P-value threshold). 
# This method does tend to miss some true positives at higher (more general) 
# levels of the GO hierarchy.
resultElim <- runTest(myGOdata.bv, algorithm="elim", statistic="fisher")
# weight01 this is the default method used by TopGO, and is a mixture of the 'elim' and 'weight' methods
resultTopgo <- runTest(myGOdata.bv, algorithm="weight01", statistic="fisher")
# when assessing a GO term, it takes into accoount the annotation of terms to the current term's parents, 
# and so reduces false positives due to the inheritance problem
resultParentchild <- runTest(myGOdata.bv, algorithm="parentchild", statistic="fisher")


# see how many results we get where weight01 gives a P-value <= 0.001:
mysummary <- summary(attributes(resultTopgo)$score <= 0.1)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001

allRes <- GenTable(myGOdata.bv, 
                   classicFisher = resultClassic, 
                   elimFisher = resultElim, 
                   topgoFisher = resultTopgo, 
                   parentchildFisher = resultParentchild, 
                   orderBy = "parentchildFisher", ranksOf = "classicFisher", topNodes = numsignif)


# write output
printGraph(myGOdata.bv, resultClassic, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultClassic", sep=""), useInfo = "all", pdfSW = TRUE)
printGraph(myGOdata.bv, resultTopgo, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultTopGo", sep=""), useInfo = "all", pdfSW = TRUE)
printGraph(myGOdata.bv, resultParentchild, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultParentchild", sep=""), useInfo = "all", pdfSW = TRUE)

write.table(allRes[,c(1,8)], file=paste(candidate_list, ".",GO_category,".GSEA_result_elimFisher.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(allRes[,c(1,7)], file=paste(candidate_list, ".",GO_category,".GSEA_result_classicFisher.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(allRes[,c(1,10)], file=paste(candidate_list, ".",GO_category,".GSEA_result_parentchild.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(allRes, file=paste(candidate_list, ".",GO_category,".GSEA_result.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.xlsx(allRes, file = paste(candidate_list, ".",GO_category,".GSEA_result.xlsx", sep=""), borders = "rows")

#### CORE MF #####
# new general run
rm(list = objects()) # clear all objects defined the previous run
#objects() # clear all objects definded the previous run

##########
# set variables
# making general script 
# Rscript myscript.R batch.csv
# and invoke these in the myscript.R
rm(myGOdata.bv,allRes,resultClassic,resultElim,resultTopgo,resultParentchild,GO_category) #clears all variables
args <- commandArgs(TRUE)
# dataset <- read.table(args[1],header=FALSE,sep=",",skip=1)

annotations=args[1]
candidate_list=args[2]
# here I will be only analyzing GO terms with at least 5 members,
# as this yield more stable results.
node_size=5

#### CORE MF #####
GO_category="MF" # molecular function
geneID2GO <- readMappings(file = annotations)  
geneUniverse <- names(geneID2GO) 

genesOfInterest.bv <- read.table(candidate_list,header=TRUE)

genesOfInterest.bv <- as.character(genesOfInterest.bv$geneID) 
geneList.bv <- factor(as.integer(geneUniverse %in% genesOfInterest.bv))
names(geneList.bv) <- geneUniverse

myGOdata.bv <- new("topGOdata", description="Candidate genes", ontology=GO_category, allGenes=geneList.bv,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = node_size)

# STATS#
# each GO term is tested independently, not taking the GO hierarchy into account
resultClassic <- runTest(myGOdata.bv, algorithm="classic", statistic="fisher")
# elim method processes the GO terms by traversing the GO hierarchy from bottom to top, 
# ie. it first assesses the most specific (bottom-most) GO terms, and proceeds later 
# to more general (higher) GO terms. When it assesses a higher (more general) GO term, 
# it discards any genes that are annotated with significantly enriched descendant 
# GO terms (considered significant using a pre-defined P-value threshold). 
# This method does tend to miss some true positives at higher (more general) 
# levels of the GO hierarchy.
resultElim <- runTest(myGOdata.bv, algorithm="elim", statistic="fisher")
# weight01 this is the default method used by TopGO, and is a mixture of the 'elim' and 'weight' methods
resultTopgo <- runTest(myGOdata.bv, algorithm="weight01", statistic="fisher")
# when assessing a GO term, it takes into accoount the annotation of terms to the current term's parents, 
# and so reduces false positives due to the inheritance problem
resultParentchild <- runTest(myGOdata.bv, algorithm="parentchild", statistic="fisher")


# see how many results we get where weight01 gives a P-value <= 0.001:
mysummary <- summary(attributes(resultTopgo)$score <= 0.1)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001

allRes <- GenTable(myGOdata.bv, 
                   classicFisher = resultClassic, 
                   elimFisher = resultElim, 
                   topgoFisher = resultTopgo, 
                   parentchildFisher = resultParentchild, 
                   orderBy = "parentchildFisher", ranksOf = "classicFisher", topNodes = numsignif)


# write output
printGraph(myGOdata.bv, resultClassic, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultClassic", sep=""), useInfo = "all", pdfSW = TRUE)
printGraph(myGOdata.bv, resultTopgo, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultTopGo", sep=""), useInfo = "all", pdfSW = TRUE)
printGraph(myGOdata.bv, resultParentchild, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultParentchild", sep=""), useInfo = "all", pdfSW = TRUE)

write.table(allRes[,c(1,8)], file=paste(candidate_list, ".",GO_category,".GSEA_result.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(allRes[,c(1,7)], file=paste(candidate_list, ".",GO_category,".GSEA_result_classicFisher.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(allRes[,c(1,10)], file=paste(candidate_list, ".",GO_category,".GSEA_result_parentchild.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(allRes, file=paste(candidate_list, ".",GO_category,".GSEA_result.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.xlsx(allRes, file = paste(candidate_list, ".",GO_category,".GSEA_result.xlsx", sep=""), borders = "rows")
