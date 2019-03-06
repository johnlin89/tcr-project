# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# BiocManager::install('rtracklayer')

library(LymphoSeq)


# Read in sample immunoSeq raw file
immunoSeq <-  read.delim("/Users/linjo/Box\ Sync/MIPs/ImmunoSeq_exported_samples/201606719-01.tsv", sep = "\t")

# Use LymphoSeq to read in all samples
TCR.list <- readImmunoSeq(path = "/Users/linjo/Box Sync/MIPS/ImmunoSeq_exported_samples")
# Filter for productive sequences
productive.TCR.nt <- productiveSeq(file.list = TCR.list,
                                   aggregate = "nucleotide", prevalence = FALSE)

alignSeq(list = productive.TCR.nt, sample = "201606719-01", type = "nucleotide", 
         method = "ClustalW", output = "consule")

lorenzCurve(samples = names(productive.TCR.nt), list = productive.TCR.nt)
topSeqsPlot(list = productive.TCR.nt, top = 10)


TCR.clonality <- clonality(file.list = productive.TCR.nt)

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