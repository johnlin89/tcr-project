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
