
library(ssh)
library(tidyverse)
library(dplyr)

# Connect to biolync server
session  <- ssh_connect("jxl2059@biolync.case.edu")

# merge with phenotypes
scp_download(session, "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink15/allWindowNoDup_flipRefPruned.fam", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
scp_download(session, "/storage/mips//MIPS_Updated.2019-02-21/data/SampleInfo/MIPs_genotyped_2019-02-15.txt", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
fam <- read.table(
  "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/allWindowNoDup_flipRefPruned.fam")
colnames(fam) <- c("FID", "IID", "PaternalID", "MaternalID", "SEX", "TCR")
eth <- read.delim(
  "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/MIPs_genotyped_2019-02-15.txt", 
  sep = "\t")
eth$IID <- paste(substring(eth$Study.ID, 1, 4), substring(eth$Study.ID, 10, 12),
                 sep = "")
eth <- select(eth, IID, Race)
eth$Race[which(eth$Race == "")] = "unknown"
eth$Race <- as.factor(eth$Race)
eth <- unique(eth)
ethFam <- merge(eth, fam, by = "IID", all.y = TRUE)
ethFam <- select(ethFam, -FID, -PaternalID,  -MaternalID)
save(ethFam, file = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/ethFam")

# Method 1: Unsupervised learning
scp_download(session, "/storage/mips//MIPS_Updated.2019-02-21/jxl2059/admixtureResults/admixture3/allWindow.2.Q", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
scp_download(session, "/storage/mips//MIPS_Updated.2019-02-21/jxl2059/admixtureResults/admixture3/crossValidation.txt", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
admixture1 <- read.table(
  "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/allWindow.2.Q")
cv <- read.table(
  "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/crossValidation.txt")
jpeg('figures/admixture1CV.jpg')
plot(cv$V4, xlab = "K", ylab = "CV Error", main = "Admixture Method 1: 10-fold Cross Validation Results", type = "l")
dev.off()
jpeg('figures/admixture1.jpg')
barplot(t(as.matrix(admixture1)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA, 
        legend.text = TRUE, main = "Admixture Method 2: Ancestry Proportions")
dev.off()

# Method 2: Supervised learning
scp_download(session, "/storage/mips//MIPS_Updated.2019-02-21/jxl2059/admixtureResults/allWindowNoDup_flipRefPruned.2.Q", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
admixture2 <- read.table(
  "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/allWindowNoDup_flipRefPruned.2.Q")
colnames(admixture2) <- c("CEU", "YRI")
# Bind with self reported races and other columns
admixture2 <- cbind(ethFam, admixture2)
# Order by percentage of race
admixture2 <- admixture2[order(admixture2$CEU), ]
# Get rid of reference populations
admixture2 <- admixture2[-grep("NA", admixture2$IID), ]
# Create column to see if results agree with self reported race
admixture2$raceTest = NA
admixture2$raceTest[which(admixture2$CEU < admixture2$YRI & 
                            admixture2$Race == "African American")] = 1
admixture2$raceTest[which(admixture2$CEU < admixture2$YRI & 
                            admixture2$Race == "European American")] = 0
jpeg('figures/admixture2.jpg')
barplot(t(as.matrix(select(admixture2, CEU, YRI))), col = rainbow(3), 
        xlab = "Individuals", ylab = "Ancestry", border = NA, axisnames = FALSE, 
        legend.text = TRUE, main = "Admixture Method 2: Ancestry Proportions")
#green is YRI
#red is CEU
dev.off()
save(admixture2, file = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/admixture2")


# Save plots to biolync server
scp_upload(session, 'figures/admixture1.jpg', to = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/figures")
scp_upload(session, 'figures/admixture2.jpg', to = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/figures")

