
library(ssh)
library(tidyverse)
library(dplyr)
# merge with phenotypes
# Connect to biolync server
session  <- ssh_connect("jxl2059@biolync.case.edu")

scp_download(session, "/storage/mips//MIPS_Updated.2019-02-21/jxl2059/admixtureResults/allWindowNoDup_flipRefPruned.2.P", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
scp_download(session, "/storage/mips//MIPS_Updated.2019-02-21/jxl2059/admixtureResults/allWindowNoDup_flipRefPruned.2.Q_se", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
scp_download(session, "/storage/mips//MIPS_Updated.2019-02-21/jxl2059/admixtureResults/allWindowNoDup_flipRefPruned.2.Q_bias", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")


# Method 1: Unsupervised learning
scp_download(session, "/storage/mips//MIPS_Updated.2019-02-21/jxl2059/admixtureResults/admixture3/allWindow.2.Q", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
scp_download(session, "/storage/mips//MIPS_Updated.2019-02-21/jxl2059/admixtureResults/admixture3/crossValidation.txt", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
admixture1 <- read.table(
  "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/allWindow.2.Q")
cv <- read.table(
  "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/crossValidation.txt")
jpeg('figures/admixture1CV.jpg')
plot(cv$V4, xlab = "K", ylab = "CV Error", main = "Admixture Method 1: 10-fold Cross Validation Results", type = "l")
line(data = cv$V4)
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
jpeg('figures/admixture2.jpg')
barplot(t(as.matrix(admixture2)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA, 
        legend.text = TRUE, main = "Admixture Method 2: Ancestry Proportions")
#green is YRI
#red is CEU
dev.off()

tbl1 = read.table("/Users/linjo/Desktop/tcr-project-desktop/snpFlipResults/windowCovar_flip.1.Q")
jpeg('figures/admixture1K1.jpg')
barplot(t(as.matrix(tbl1)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA)
dev.off()
tbl2 = read.table("/Users/linjo/Desktop/tcr-project-desktop/snpFlipResults/windowCovar_flip.2.Q")
jpeg('figures/admixture1K2.jpg')
barplot(t(as.matrix(tbl2)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA)
dev.off()
tbl3 = read.table("/Users/linjo/Desktop/tcr-project-desktop/snpFlipResults/windowCovar_flip.3.Q")
jpeg('figures/admixture1K3.jpg')
barplot(t(as.matrix(tbl3)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA)
dev.off()
tbl4 = read.table("/Users/linjo/Desktop/tcr-project-desktop/snpFlipResults/windowCovar_flip.4.Q")
jpeg('figures/admixture1K4.jpg')
barplot(t(as.matrix(tbl4)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA)
dev.off()
tbl5 = read.table("/Users/linjo/Desktop/tcr-project-desktop/snpFlipResults/windowCovar_flip.5.Q")
jpeg('figures/admixture1K5.jpg')
barplot(t(as.matrix(tbl5)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA)
dev.off()


### all samples
tbl6 = read.table("/Users/linjo/Desktop/tcr-project-desktop/snpFlipResults/allWindow_flip.5.Q")
jpeg('figures/allAdmixture1K1.jpg')
barplot(t(as.matrix(tbl6)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA)
dev.off()
