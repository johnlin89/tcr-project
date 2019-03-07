# Copyright: CC-BY-SA 4.0
# Author Notes:  Please note this work is in conjunction with Dr. Dana Crawford and Dr. William Bush at Case Western Reserve Universithy and Cleveland Institute of Computational Biology. The data used below is not my own.
# File Description: The code below is intended to visualize the immunoSeq data

# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# BiocManager::install('rtracklayer')

library(LymphoSeq)

# TODO: Need to store on the biolync server with RCurl
# Use LymphoSeq to read in all samples
TCR.list <- readImmunoSeq(path = 
            "/Users/linjo/Box Sync/MIPS/ImmunoSeq_exported_samples")

# Filter for productive sequences
productive.TCR.nt <- productiveSeq(file.list = TCR.list,
                                   aggregate = "nucleotide", prevalence = FALSE)

# Generate and save visualization plots
# TODO: use Rcurl to read this in from server
jpeg('lorenz1.jpg')
lorenzCurve(samples = names(productive.TCR.nt), list = productive.TCR.nt)
dev.off()
jpeg('topSeq1.jpb')
topSeqsPlot(list = productive.TCR.nt, top = 10)
dev.off()

### Experimental Code

alignSeq(list = productive.TCR.nt, sample = "201606719-01", type = "nucleotide", 
         method = "ClustalW", output = "consule")




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
