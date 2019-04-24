library(tidyverse)
library(dplyr)
# merge with phenotypes
# Connect to biolync server
session  <- ssh_connect("jxl2059@biolync.case.edu")

eigenv <- read.table('/Users/linjo/Desktop/tcr-project-desktop/pca/pcaResult.eigenvec')


scp_download(session, "/storage/mips//MIPS_Updated.2019-02-21/data/SampleInfo/MIPs_genotyped_2019-02-15.txt", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
eth <- read.delim(
  "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/MIPs_genotyped_2019-02-15.txt", 
  sep = "\t")
eth$IID <- paste(substring(eth$Study.ID, 1, 4), substring(eth$Study.ID, 10, 12),
                                                        sep = "")

eth <- select(eth, IID, Race)
colnames(eigenv)[2] <- "IID"

test <- merge(eth, eigenv, by = "IID")
test$color = as.integer(test[,2])
pairs(test[,3:7], col=test$color, pch=20, lower.panel=NULL)
