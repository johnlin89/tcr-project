# Copyright: CC-BY-SA 4.0
# Author Notes:  Please note this work is in conjunction with Dr. Dana Crawford and Dr. William Bush at Case Western Reserve Universithy and Cleveland Institute of Computational Biology. The data used below is not my own.
# File Description: The code below is intended to clean and format the phenotype data in order to combine with the germline genotype data in plink.

require(tidyverse)
require(RCurl)

# Create phenotype file tcrEmrPheno.txt to input into plink
# Read in all Phenotypes
# TODO: Need to store on the biolync server with RCurl
tcrEmrPheno <- read.delim(
  "/Users/linjo/Box\ Sync/MIPs/Immunoseq_EHR_2019-02-15_for_BOX.txt", 
  sep = "\t", strip.white = TRUE, stringsAsFactors = FALSE, na.strings = "")
# Create FIID and IID columns
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
# Output for plink
# TODO: Need to store on the biolync server with RCurl
write_delim(tcrEmrPheno,
            path = "/Users/linjo/Box\ Sync/MIPs/tcrEmrPheno.txt",
            delim = "\t", col_names = TRUE, quote_escape = FALSE, na = "missing")
