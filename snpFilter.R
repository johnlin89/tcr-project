# Copyright statement comment
# Author comment
# File description comment, including purpose of program, inputs, and outputs
# source() and library() statements

# install.packages("KRIS")
# install.packages("pheatmap")
# install.packages("RCurl")
# install.packages("qqman")
library(tidyverse)
library(rtracklayer)
library(KRIS)
library(pheatmap)
library(LymphoSeq)
library(RCurl)
# Manhattan plots
library(qqman)

# Look at EMR Data
emrAll <- read.delim("/Users/linjo/Box\ Sync/MIPs/Biorepository-Data.V2_2017-06-01_for_BOX.txt",
                     sep = "\t")
emr15 <- read.delim("/Users/linjo/Box\ Sync/MIPs/Immunoseq_EHR_2019-02-15_for_BOX.txt",
                    sep = "\t")

# Create phenotype file tcrEmrPheno.txt to input into plink
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
            path = "/Users/linjo/Box\ Sync/MIPs/tcrEmrPheno.txt",
            delim = "\t", col_names = TRUE, quote_escape = FALSE)

# Define snp window region
# http://www.imgt.org/IMGTrepertoire/LocusGenes/chromosomes/human/Hu_IGTRloci.html
tcrbChr <- 7
# start and end in Mb
tcrbChrStart <- 141.2
tcrbChrEnd <- 142

# Read in plink files
plinkRead <- read.bed("/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/jxl2059/plinkResults/linearTest/MIPS_SexCorrected_Pheno_tcrb_linear.bed", "/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/jxl2059/plinkResults/linearTest/MIPS_SexCorrected_Pheno_tcrb_linear.bim", "/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/jxl2059/plinkResults/linearTest/MIPS_SexCorrected_Pheno_tcrb_linear.fam", only.snp = FALSE)
#http://zzz.bwh.harvard.edu/plink/binary.shtml
#2 is homozygote A2A2
#1 is hetero A1A2
#0 is homozygote?? A1A1
#NA is missing

plinkPed <- plinkRead$snp
plinkMap <- plinkRead$snp.info
plinkFam <- plinkRead$ind.info

plinkPed <- as.data.frame(plinkPed)

# Add SNP column names to plinkPed, is order same??
colnames(plinkPed) <- plinkMap$ID
colnames(plinkPed)[grep(",", colnames(plinkPed))] <- gsub(",","|",colnames(plinkPed)[grep(",", colnames(plinkPed))])

# Filter plinkPed for 15 patients with TCR metrics, verify order is ok??
# Also sdd in phenotypes
plinkPed$FID <- plinkFam$FamID
plinkPed$IID <- plinkFam$IndID
plinkPed$phenotype <- plinkFam$phenotype
plinkPed <- plinkPed[which(plinkPed$IID %in% tcrEmrPheno$IID),]

# Filter plinkPed for relevant snps?? how should we filter snps??
plinkPcMeans <- read.table("/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/jxl2059/plinkResults/linearTest/MIPS_SexCorrected_Pheno_tcrb_linear.qassoc.means", header = TRUE)
# only grab where counts are 5 for each category
snpFiltered <- plinkPcMeans$SNP[which(plinkPcMeans$VALUE == "COUNTS" & plinkPcMeans$G11 == 5 & plinkPcMeans$G12 == 5 & plinkPcMeans$G11 == 5)]
plinkPed <- plinkPed[,which(colnames(plinkPed) %in% c(as.character(snpFiltered), "FID", "IID", "phenotype"))]

#colnames(plinkPed)[(!(colnames(plinkPed) %in% c("FID", "IID", "phenotype")))]

for (i in 1:length(snpFiltered)) {
  if (i == 1) {
    snpFilteredList <- as.character(snpFiltered[i])
  } else {
    snpFilteredList <- paste(snpFilteredList, as.character(snpFiltered[i]), sep = ", ")
  }
}

by_snp <- group_by_(plinkPed, .dots = names(plinkPed)[-grep("FID|IID|phenotype", names(plinkPed))])
by_snp_means <- summarise(by_snp, meanTcr = mean(phenotype, na.rm = TRUE))

#####playing with the below


plinkLinear <- read.table("/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/jxl2059/plinkResults/linearTest/MIPS_SexCorrected_Pheno_tcrb_linear.assoc.linear", header = TRUE)
# Filter for snps with 15 people
plinkLinear <-  plinkLinear[which(plinkLinear$NMISS == 15),]
plinkLinear <- plinkLinear[which(plinkLinear$P <= 0.05),]

plinkPcMeans <- plinkPcMeans[which(plinkPcMeans$SNP %in% snpFiltered),]

tes <- reshape(plinkPcMeans, v.names = c("G11", "G12", "G22"), idvar = "SNP", timevar = "VALUE", direction = "wide")

test <- spread(plinkPcMeans, key = "VALUE", value = c("G11", "G12", "G22"))

manhattan(plinkLinear, ylim = c(0, 10))
qq(plinkLinear$P)
