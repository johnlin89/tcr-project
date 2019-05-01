library(tidyverse)
library(dplyr)
library(ggplot2)
library(GGally)
# merge with phenotypes
# Connect to biolync server
session  <- ssh_connect("jxl2059@biolync.case.edu")

eigenv <- read.table('/Users/linjo/Desktop/tcr-project-desktop/pca/pcaResult.eigenvec', header = FALSE, skip = 0, sep = " ")

scp_download(session, "/storage/mips//MIPS_Updated.2019-02-21/data/SampleInfo/MIPs_genotyped_2019-02-15.txt", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
eth <- read.delim(
  "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/MIPs_genotyped_2019-02-15.txt", 
  sep = "\t")
eth$IID <- paste(substring(eth$Study.ID, 1, 4), substring(eth$Study.ID, 10, 12),
                 sep = "")

eth <- select(eth, IID, Race)
eth$Race[which(eth$Race == "")] = "unknown"


eth$Race <- as.factor(eth$Race)
eth <- unique(eth)
#Read in the eigenvectors
eigenv <- data.frame(eigenv)
colnames(eigenv)[2] <- "IID"
eigenMerge <- merge(eth, eigenv, by = "IID")
rownames(eigenMerge) <- eigenMerge$IID

eigenMerge <- select(eigenMerge, -IID, -V1)

eigenMerge$color = as.integer(eigenMerge$Race)

colnames(eigenMerge)[2:21] = paste("PC", 1:20, sep = "")

pairs(eigenMerge[, 2:6], col = eigenMerge$color, pch = 20, lower.panel = NULL)
# need to switch this on and off
par(xpd = TRUE)
legend("bottomleft", fill = unique(eigenMerge$Race), legend = c(levels(eigenMerge$Race)), cex = 0.50)

eigenMerge2 <- select(eigenMerge, -Race)
proportionvariances <- ((apply(eigenMerge2, 1, sd)^2) / (sum(apply(eigenMerge2, 1, sd)^2)))*100



# plot variance
eigenval <- read.table('/Users/linjo/Desktop/tcr-project-desktop/pca/pcaResult.eigenval', header = FALSE, skip = 0, sep = " ")
eigenval$var <- eigenval$V1/sum(eigenval$V1)
eigenval$var2 <- cumsum(eigenval$var)
jpeg('figures/pcaVariance.jpg')
ggplot() + 
  geom_point(data = data.frame(eigenval), mapping = aes(x = 1:20, y = var), col = "red") +
  geom_line(data = data.frame(eigenval), mapping = aes(x = 1:20, y = var), col = "blue") +
  xlab("Principle Components") + ylab("% Variance") + 
  ggtitle("Variance Explained by Individual Principle Components")
dev.off()