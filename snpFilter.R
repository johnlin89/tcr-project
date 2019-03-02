# Copyright statement comment
# Author comment
# File description comment, including purpose of program, inputs, and outputs
# source() and library() statements

# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# BiocManager::install('rtracklayer')

#install.packages("KRIS")

#install.packages("pheatmap")
library(tidyverse)
library(rtracklayer)
library(KRIS)
library(pheatmap)
library(LymphoSeq)

# Read in plink files
plinkRead <- read.bed("/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected.bed", "/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected.bim", "/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected.fam", only.snp = FALSE)

# Read in all Phenotypes
tcrEmrPheno <- read.delim(
  "/Users/linjo/Box\ Sync/MIPs/Immunoseq_EHR_2019-02-15_for_BOX.txt", 
  sep = "\t", strip.white = TRUE)
# Create FIID and IID columns
tcrEmrPheno$FID <- paste(substring(tcrEmrPheno$MIPs.ID, 0, 4), 
                         substring(tcrEmrPheno$MIPs.ID, 10, 13), sep = "")
tcrEmrPheno$IID <- paste(substring(tcrEmrPheno$MIPs.ID, 0, 4), 
                         substring(tcrEmrPheno$MIPs.ID, 10, 13), sep = "")
tcrEmrPheno <- tcrEmrPheno[, c(47:48, 1:46)]
# Replace unwanted characters
tcrEmrPheno <- data.frame(lapply(tcrEmrPheno, function(x) {
               gsub(" ", "_", x)
               }))
tcrEmrPheno <- data.frame(lapply(tcrEmrPheno, function(x) {
               gsub(",", "", x)
               }))
tcrEmrPheno <- data.frame(lapply(tcrEmrPheno, function(x) {
               gsub("%", "", x)
               }))
# Output for plink
write_delim(tcrEmrPheno,
            path = "/Users/linjo/Box\ Sync/MIPs/20190301_tcrEmrPheno.txt",
            delim = "\t", col_names = TRUE, quote_escape = FALSE)

# Define snp window region
# http://www.imgt.org/IMGTrepertoire/LocusGenes/chromosomes/human/Hu_IGTRloci.html
tcrbChr <- 7
tcrbChrStart <- 141.2
tcrbChrEnd <- 142






#http://zzz.bwh.harvard.edu/plink/binary.shtml

#2 is homozygote A2A2
#1 is hetero A1A2
#0 is homozygote?? A1A1
#NA is missing
emr <- read.delim("/Users/linjo/Box\ Sync/MIPs/Biorepository-Data.V2_2017-06-01_for_BOX.txt",
                  sep = "\t")
emr2 <- read.delim("/Users/linjo/Box\ Sync/MIPs/Immunoseq_EHR_2019-02-15_for_BOX.txt",
                   sep = "\t")
immunoSeq <-  read.delim("/Users/linjo/Box\ Sync/MIPs/ImmunoSeq_exported_samples/201606719-01.tsv",
                         sep = "\t")







productive.TCR.nt <- productiveSeq(file.list = TCR.list,
                                   aggregate = "nucleotide", prevalence = FALSE)

TCR.clonality <- clonality(file.list = productive.TCR.nt)

TCR.list <- readImmunoSeq(path = "/Users/linjo/Box Sync/MIPS/ImmunoSeq_exported_samples")
vGenes <- geneFreq(productive.nt = productive.TCR.nt, locus = "V", family = TRUE)
library(reshape)
vGenes <- reshape::cast(vGenes, familyName ~ samples, value = "frequencyGene", sum)
rownames(vGenes) = as.character(vGenes$familyName)
vGenes$familyName = NULL
RedBlue <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(256)
pheatmap::pheatmap(vGenes, color = RedBlue, scale = "row")


clonalRelatedness(list = TCR.list, editDistance = 10)


phyloTree(list = TCR.list, sample = "201606719-01", type = "nucleotide", 
          layout = "rectangular")
