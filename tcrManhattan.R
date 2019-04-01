# Copyright: CC-BY-SA 4.0
# Author Notes:  Please note this work is in conjunction with Dr. Dana Crawford and Dr. William Bush at Case Western Reserve Universithy and Cleveland Institute of Computational Biology. The data used below is not my own.
# File Description: The code below is intended to work with the output from phenoClean.R and the sex-corrected genotype data on the biolync.case.edu server. The goal is to create a Manhattan plot delineating significant SNPs.

require(KRIS)
require(RCurl)
require(rtracklayer)
require(pheatmap)
require(qqman)
require(ssh)

# Connect to biolync server
session <- ssh_connect("jxl2059@biolync.case.edu")

# Demonstration of single SNP test

# Perform linear regression for 1 SNP
snp <- "rs11768792"
plinkPedSingle <- plinkPed[, c(snp,"phenotype")]
lm.fit <- lm(phenotype ~ rs11768792, data = plinkPedSingle)
jpeg("linearRegressionExample.jpg")
ggplot(data = plinkPedSingle) + 
  geom_point(mapping = aes(x = rs11768792, y = phenotype), color = "blue") +
  geom_line(aes(x = rs11768792, y = predict(lm.fit)), color = "red") + 
  xlab(paste("snp", snp)) + ylab("TCR Productive Clonality") +
  ggtitle("SNP Genotype vs TCR Productive Clonality")
dev.off()
summary(lm.fit)
plinkLinear[which(plinkLinear$SNP == "rs11768792"), ]


# Zoomed manhattan with 50,000 bp window on each side
# Read in linear regression results
plinkLinear <- scp_download(session, "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink5/window.assoc.linear", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
plinkLinear <- read.table("/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/window.assoc.linear", header = TRUE)
# Remove where snps with P values of NA
plinkLinear <- na.omit(plinkLinear, col = "P")
# Adjust p-values
# p.adjust(plinkLinear$P, method = "bonferroni")
# Generate and save Manhattan plot
jpeg('figures/manhattanZoomed.jpg')
# coordinates for expanded window are 141948851, 142560972
manhattan(plinkLinear, xlim = c(141900000, 142660000), 
          suggestiveline = -log10(0.05), annotatePval = 0.05,
          highlight = as.character(plinkLinear$SNP[which(plinkLinear$P < 0.05)]))
dev.off()

# Print significant snps
print("These SNPs appear to be significant...")
plinkLinear[which(plinkLinear$P < 0.05),]

# Generate and save qq plot
jpeg('figures/qqplotZoomed.jpg')
qq(plinkLinear$P)
dev.off()
# Save plots to biolync server
scp_upload(session, 'figures/manhattanZoomed.jpg', to = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/figures")
scp_upload(session, 'figures/qqplotZoomed.jpg', to = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/figures")

# All of chr 7 manhattan
plinkLinearAll <- scp_download(session, "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink6/allChr7.assoc.linear", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
plinkLinearAll <- read.table("/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/allChr7.assoc.linear", header = TRUE)
# Remove where snps with P values of NA
plinkLinearAll <- na.omit(plinkLinearAll, col = "P")
#Adjust p-values
# p.adjust(plinkLinear$P, method = "bonferroni")
# Generate and save Manhattan plot
jpeg('figures/manhattanAll.jpg')
manhattan(plinkLinearAll, suggestiveline = -log10(0.05))
dev.off()
# Generate and save qq plot
jpeg('figures/qqplotAll.jpg')
qq(plinkLinearAll$P)
dev.off()
# Save plots to biolync server
scp_upload(session, 'figures/manhattanAll.jpg', to = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/figures")
scp_upload(session, 'figures/qqplotAll.jpg', to = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/figures")

# Adding in covariates for age and gender
plinkLinearCovar <- scp_download(session, "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink7/windowCovar.assoc.linear", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
plinkLinearCovar <- read.table("/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/windowCovar.assoc.linear", header = TRUE)
# Remove where snps with P values of NA
plinkLinearCovar <- na.omit(plinkLinearCovar, col = "P")
# Remove P-values not associated with SNPs
plinkLinearCovar <- plinkLinearCovar[which(plinkLinearCovar$TEST == "ADD"),]
# Adjust p-values
# p.adjust(plinkLinear$P, method = "bonferroni")
# Generate and save Manhattan plot
jpeg('figures/manhattanCovar.jpg')
manhattan(plinkLinearCovar, xlim = c(141948851, 142560972), 
          suggestiveline = -log10(0.05), annotatePval = 0.05,
          highlight = as.character(plinkLinearCovar$SNP[which(plinkLinearCovar$P < 0.05)]))
dev.off()

# Print significant snps
print("These SNPs appear to be significant...")
plinkLinearCovar[which(plinkLinearCovar$P < 0.05),]

# Generate and save qq plot
jpeg('figures/qqplotCovar.jpg')
qq(plinkLinearCovar$P)
dev.off()
# Save plots to biolync server
scp_upload(session, 'figures/manhattanCovar.jpg', to = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/figures")
scp_upload(session, 'figures/qqplotCovar.jpg', to = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/figures")


# Disconnect session
ssh_disconnect(session)

### Experimental Code

# Read in plink files
# http://zzz.bwh.harvard.edu/plink/binary.shtml
# TODO: Need to store on the biolync server with RCurl
plinkRead <- read.bed("/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/jxl2059/plinkResults/linearTest/MIPS_SexCorrected_Pheno_tcrb_linear.bed", "/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/jxl2059/plinkResults/linearTest/MIPS_SexCorrected_Pheno_tcrb_linear.bim", "/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/jxl2059/plinkResults/linearTest/MIPS_SexCorrected_Pheno_tcrb_linear.fam", 
                      only.snp = FALSE)

# Creating separate data frames from bed object
plinkPed <- plinkRead$snp
plinkPed <- as.data.frame(plinkPed)
plinkMap <- plinkRead$snp.info
plinkFam <- plinkRead$ind.info

# Add SNP column names to plinkPed
# TODO: Verify order of snps from plinkMap can be transposed
colnames(plinkPed) <- plinkMap$ID
# Substitute commas for pipes in column names of snps
colnames(plinkPed)[grep(",", colnames(plinkPed))] <- 
  gsub(",","|",colnames(plinkPed)[grep(",", colnames(plinkPed))])

# Filter plinkPed for 15 patients with TCR metrics
# TODO: Verify we can just bind the plinkFam columns - the order is the same?
plinkPed$FID <- plinkFam$FamID
plinkPed$IID <- plinkFam$IndID
# Also, add in phenotypes to plinkPed
# TODO: Verify - the order is the same?
plinkPed$phenotype <- plinkFam$phenotype
# Filter plinkPed for 15 patients with TCR metrics
plinkPed <- plinkPed[which(plinkPed$IID %in% tcrEmrPheno$IID),]

# Read in productive clonality means associated with each haplotype
plinkPcMeans <- read.table("/Users/linjo/Desktop/tcr-project-desktop/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink1.qassoc.means", 
                           header = TRUE)
# Filter plinkPed for relevant snps
# Only grab where counts are 5 for each category
snpFiltered <- plinkPcMeans$SNP[which(plinkPcMeans$VALUE == "COUNTS" & 
                                        plinkPcMeans$G11 == 5 & 
                                        plinkPcMeans$G12 == 5 & 
                                        plinkPcMeans$G11 == 5)]
plinkPed <- plinkPed[,which(colnames(plinkPed) %in% 
                              c(as.character(snpFiltered), 
                                "FID", "IID", "phenotype"))]
plinkPed <- select(plinkPed, FID, IID, phenotype, everything())
for (i in 1:length(snpFiltered)) {
  if (i == 1) {
    snpFilteredList <- as.character(snpFiltered[i])
  } else {
    snpFilteredList <- paste(snpFilteredList, as.character(snpFiltered[i]), sep = ", ")
  }
}
# Manually determine means
by_snp <- group_by_(plinkPed, .dots = names(plinkPed)[-grep("FID|IID|phenotype", names(plinkPed))])
by_snp_means <- summarise(by_snp, meanTcr = mean(phenotype, na.rm = TRUE))

# Demonstration of single SNP test
snp <- "rs11768792"
plinkPedSingle <- plinkPed[, c(snp,"phenotype")]
lm.fit <- lm(phenotype ~ rs6945601, data = plinkPedSingle)
jpeg("linearRegressionExample.jpg")
ggplot(data = plinkPedSingle) + 
  geom_point(mapping = aes(x = rs6945601, y = phenotype), color = "blue") +
  geom_line(aes(x = rs6945601, y = predict(lm.fit)), color = "red") + 
  xlab(paste("snp", snp)) + ylab("TCR Productive Clonality") +
  ggtitle("SNP Genotype vs TCR Productive Clonality")
dev.off()
summary(lm.fit)
plinkLinear[which(plinkLinear$SNP == "rs6945601"), ]


