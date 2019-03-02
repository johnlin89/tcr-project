# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# BiocManager::install('rtracklayer')

#install.packages("KRIS")

install.packages("pheatmap")
library(rtracklayer)
library(KRIS)
library(pheatmap)
library(LymphoSeq)


tcrBStart <- 142299011




# test2 <- read.bed("/Users/linjo/GoogleDrive/CaseWesternUniversity/ckd-project/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected.bed", "/Users/linjo/GoogleDrive/CaseWesternUniversity/ckd-project/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected.bim", "/Users/linjo/GoogleDrive/CaseWesternUniversity/ckd-project/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected.fam", only.snp = FALSE)

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
