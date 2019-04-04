# Copyright: CC-BY-SA 4.0
# Author Notes:  Please note this work is in conjunction with Dr. Dana Crawford and Dr. William Bush at Case Western Reserve Universithy and Cleveland Institute of Computational Biology. The data used below is not my own.
# File Description: The code below is intended to clean and format the phenotype data in order to combine with the germline genotype data in plink.

require(tidyverse)
require(RCurl)
require(ssh)

# Connect to biolync server
session <- ssh_connect("jxl2059@biolync.case.edu")

# Create phenotype file tcrEmrPheno.txt to input into plink
# Read in all Phenotypes
scp_download(session, "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/phenotypes/Immunoseq_EHR_2019-02-15_for_BOX.txt", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
tcrEmrPheno <- read.delim(
  "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/Immunoseq_EHR_2019-02-15_for_BOX.txt", 
  sep = "\t", strip.white = TRUE, stringsAsFactors = FALSE, na.strings = "")
# Create FID and IID columns
tcrEmrPheno$FID <- paste(substring(tcrEmrPheno$MIPs.ID, 0, 4), 
                         substring(tcrEmrPheno$MIPs.ID, 10, 13), sep = "")
tcrEmrPheno$IID <- paste(substring(tcrEmrPheno$MIPs.ID, 0, 4), 
                         substring(tcrEmrPheno$MIPs.ID, 10, 13), sep = "")
# Rearrange so that FID, IID, and prod clon are first
tcrEmrPheno <- select(tcrEmrPheno, FID, IID, Productive.Clonality, everything())
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
# Rename Gender to SEX and Age to AGE
colnames(tcrEmrPheno)[which(colnames(tcrEmrPheno) == "Gender")] = "SEX"
colnames(tcrEmrPheno)[which(colnames(tcrEmrPheno) == "Age")] = "AGE"
# Record SEX Coded 1/2/0 for M/F/missing per PLINK
tcrEmrPheno$SEX = as.character(tcrEmrPheno$SEX)
tcrEmrPheno$SEX[which(tcrEmrPheno$SEX == "M")] = 1
tcrEmrPheno$SEX[which(tcrEmrPheno$SEX == "F")] = 2
# Filter for columns of interest
# tcrEmrPheno <- select(tcrEmrPheno, FID, IID, Productive.Clonality, SEX, AGE)

# Data Visualization and EDA
# Create density plots and histograms
jpeg('figures/histTcr.jpg')
hist(as.numeric(as.character(tcrEmrPheno$Productive.Clonality)), 
     main = "Histogram Distribution of Productive Clonality (n = 15)",
     xlab = "Productive Clonality",
     ylab = "Number of individuals")
dev.off()
jpeg('figures/densTcr.jpg')
plot(density(as.numeric(as.character(tcrEmrPheno$Productive.Clonality))),
     main = "Density Plot of Productive Clonality (n = 15)",
     xlab = "Productive Clonality",
     ylab = "Number of individuals")
dev.off()

# Output for plink and store on biolync
scp_upload(session, 'figures/densTcr.jpg', to = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/figures")
scp_upload(session, 'figures/histTcr.jpg', to = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/figures")
write_delim(tcrEmrPheno,
            path = "/Users/linjo/Desktop/tcr-project-desktop/tcrEmrPheno.txt",
            delim = "\t", col_names = TRUE, quote_escape = FALSE, na = "missing")
scp_upload(session, files = "/Users/linjo/Desktop/tcr-project-desktop/tcrEmrPheno.txt", to = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/phenotypes")

# Disconnect session
ssh_disconnect(session)

