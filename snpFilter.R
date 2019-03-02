# Copyright statement comment
# Author comment
# File description comment, including purpose of program, inputs, and outputs
# source() and library() statements

#install.packages("KRIS")
#install.packages("pheatmap")
#install.packages("RCurl")
library(tidyverse)
library(rtracklayer)
library(KRIS)
library(pheatmap)
library(LymphoSeq)
library(RCurl)

# Read in plink files
plinkRead <- read.bed("/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected.bed", "/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected.bim", "/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected.fam", only.snp = FALSE)
#http://zzz.bwh.harvard.edu/plink/binary.shtml
#2 is homozygote A2A2
#1 is hetero A1A2
#0 is homozygote?? A1A1
#NA is missing


# Look at EMR Data
emrAll <- read.delim("/Users/linjo/Box\ Sync/MIPs/Biorepository-Data.V2_2017-06-01_for_BOX.txt",
                     sep = "\t")
emr15 <- read.delim("/Users/linjo/Box\ Sync/MIPs/Immunoseq_EHR_2019-02-15_for_BOX.txt",
                    sep = "\t")

# Create phenotype file 20190301_tcrEmrPheno.txt to input into plink
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
# start and end in Mb
tcrbChrStart <- 141.2
tcrbChrEnd <- 142

plinkLinear <- read.table("/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/jxl2059/plinkResults/linearTest/MIPS_SexCorrected_Pheno_tcrb_linear.assoc.linear", header = TRUE)
# Filter for snps with 15 people
plinkLinear <-  plinkLinear[which(plinkLinear$NMISS == 15),]
plinkLinear <- plinkLinear[which(plinkLinear$P <= 0.05),]




